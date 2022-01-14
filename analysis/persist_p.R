library(data.table)
library(limma)
library(errors)
options(stringsAsFactors=F)

#Read in data
constructs = data.frame(fread("./constructs.tsv"))
rownames(constructs) = constructs$SequenceID
samples = data.frame(fread("./samples.csv",header=T))
rownames(samples) = samples$index
dat = data.frame(fread("./joined.txt",header=T))
counts = as.matrix(dat[, -1])
rownames(counts) = dat$id
colnames(counts) = unlist(lapply(strsplit(colnames(counts),"\\."), "[[", 4))
spikeins = c("120002B1", "120010B1", "220023B1")

#In-cell stability data processing
ics_samples = subset(samples, 
                     (Experiment=="stability")&
                     (Spikein==1)&
                     (Group=="In-cell")&
                     (Type=="Long")
                    )$index
ics_labels = samples[ics_samples, "Time"]
ics = as.data.frame(log2(counts[c(subset(constructs, longamp)$SequenceID, spikeins) ,ics_samples]))
colnames(ics) = ics_labels

#Scaling with spike-ins
ics_scaled=ics
for (i in 1:ncol(ics))
{
    spikeinfit = lm(ics[spikeins, c(1, i)])
    ics_scaled[, i]=predict.lm(spikeinfit, newdata=ics[, c(1, i)])
}

#Calculate degradation rates; uses lmFit function from limma
ICS_USED = c(1, 7, 12) #Drop later times
ics_design = model.matrix(~ics_labels[ics_labels%in%ICS_USED])
colnames(ics_design)[2] = "time"
v_ics = voom(2^ics_scaled[!rownames(ics_scaled)%in%spikeins, ics_labels%in%ICS_USED],
             design=ics_design, lib.size=rep(1, ncol(ics_scaled[, ICS_USED])))
fit_ics = lmFit(v_ics, ics_design)
coef_ics = fit_ics$coefficients
se_ics = fit_ics$stdev.unscaled*fit_ics$sigma #Unmod. standard error from reg.

degrate_ics = -coef_ics[, 2]
error_ics = se_ics[, 2]
int_ics = coef_ics[, 1]
halflife_ics = log(2) / degrate_ics

#Polysome data processing
FRACTION_START = 2
FRACTION_END = 16
WEIGHT_FRACTION_END = 9
FRACTION_WEIGHTS = c(0,0,0,1,1.6,2.7,3.8,4,5.5,7,9,11.5,14,18,24,30)
MONOSOME_FRACTION = "4"
PRE_FRACTION = "2"
POLY_RANGE = as.character(7:9)

poly_samples = subset(samples,(Experiment=="polysomes")&(Spikein==1))$index
poly_labels = samples[poly_samples,"Fraction"]
poly = as.data.frame(log2(counts[, poly_samples]))
colnames(poly) = poly_labels

#Scaling with spikeins
poly_scaled=poly
for (i in 1:ncol(poly))
{
    spikeinfit = lm(poly[spikeins, c(1, i)])
    poly_scaled[,i] = predict.lm(spikeinfit, newdata=poly[, c(1, i)])
}

poly_fraction = 2^poly_scaled[, FRACTION_START:FRACTION_END] / rowSums(2^poly_scaled[, FRACTION_START:FRACTION_END])
poly_fraction_weighted = t(t(poly_fraction)*FRACTION_WEIGHTS[FRACTION_START:FRACTION_END])
poly_ribonum = rowSums(poly_fraction_weighted[, 1:(WEIGHT_FRACTION_END-FRACTION_START+1)])

#Ribonum bootstrapping
poly_nsample = 10000
poly_sampled_ribonum_all = matrix(nrow=length(poly_ribonum),ncol=poly_nsample, NA)
for(ii in 1:poly_nsample)
{
    poly_scaled_sampled=poly
    for (i in 1:ncol(poly))
    {
        spikeinfit = lm(poly[spikeins, c(1, i)])
        sampled_fractions = sample(1:ncol(poly), 1)
        sampled_data = poly[, c(1,sampled_fractions)]
        colnames(sampled_data)[2] = colnames(spikeinfit$model)[2]
        poly_scaled_sampled[,i] = predict.lm(spikeinfit, newdata=sampled_data)
    }
    
    poly_sampled_fraction = 2^poly_scaled_sampled[, FRACTION_START:FRACTION_END] / rowSums(2^poly_scaled[, FRACTION_START:FRACTION_END])
    poly_sampled_fraction_weighted = t(t(poly_sampled_fraction)*FRACTION_WEIGHTS[FRACTION_START:FRACTION_END])
    poly_sampled_ribonum = rowSums(poly_sampled_fraction_weighted[, 1:(WEIGHT_FRACTION_END-FRACTION_START+1)])
    poly_sampled_ribonum_all[, ii] = poly_sampled_ribonum
}
poly_ribonum_error = apply(poly_sampled_ribonum_all, 1, sd)
names(poly_ribonum_error) = names(poly_ribonum)

#Fraction ratios
ratio_mono_pre = poly_fraction[,MONOSOME_FRACTION] / poly_fraction[,PRE_FRACTION]
names(ratio_mono_pre) = rownames(poly_fraction)
ratio_poly_mono = rowSums(poly_fraction[,POLY_RANGE]) / poly_fraction[,MONOSOME_FRACTION]
names(ratio_poly_mono) = rownames(poly_fraction)
ratio_poly_pre = rowSums(poly_fraction[,POLY_RANGE]) / poly_fraction[,PRE_FRACTION]
names(ratio_poly_pre) = rownames(poly_fraction)

used_ids = intersect(names(degrate_ics), names(poly_ribonum))
out = data.frame(cbind(do.call(cbind, lapply(list(degrate_ics, error_ics, poly_ribonum, poly_ribonum_error), "[", used_ids))))
colnames(out) = c("degrate", "degrate_err", "ribonum", "ribonum_err")
rownames(out) = used_ids

#Expected protein level
kdeg_prot = log(2) / Inf #>6h From Promega
prot_exp = list()
times = seq(0, 64, by=0.5) #Time steps
for(id in rownames(out))
{
    inv_mrna_length = 1 / constructs[id, "effective_mw"] #effective MW with +150A's and AUGC content
    num_ribosome = poly_ribonum[id]
    errors(num_ribosome) = poly_ribonum_error[id]
    kdeg_mrna = degrate_ics[id]
    errors(kdeg_mrna) = error_ics[id]
    prot_exp[[id]] = num_ribosome * inv_mrna_length / (kdeg_prot - kdeg_mrna) * (exp(-kdeg_mrna * times) - exp(-kdeg_prot * times))
}  

prot_exp_mat = drop_errors(t(do.call(cbind, prot_exp)))
prot_exp_mat_err = t(do.call(cbind, lapply(prot_exp, errors)))
colnames(prot_exp_mat) = times
colnames(prot_exp_mat_err) = times

#Dump output
write.table(poly_fraction, "poly_fractions.tsv", row.names=T, col.names=T, quote=F, sep="\t")
write.table(out, "degrate_ribonum.tsv", row.names=T, col.names=T, quote=F, sep="\t")
write.table(prot_exp_mat, "expression_pred.tsv", row.names=T, col.names=T, quote=F, sep="\t")
write.table(prot_exp_mat, "expression_pred_err.tsv", row.names=T, col.names=T, quote=F, sep="\t")
