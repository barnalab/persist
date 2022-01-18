%% Read in  Gun's processed counts -- UPDATED with new barcodes
x = readtable( 'GunAnalysis_ivstability/long_spikein_invitrostability.txt' );
dx = table2cell( x );
barcodes_ivs = dx(:,1);
log2d = cell2mat(dx(:,2:end));
d = 2.^log2d;
times = [0,0.5,1,2,3,4,5,6,16,24]*60; % in hours.
d0 = d(:,1);
d_norm = d./repmat( d0, 1, size(d,2));
semilogy( times/60, d_norm','o-' );
xlabel( 'Time (hours)' ); ylabel( 'Frac. full-length RT-able' );
set(gcf, 'PaperPositionMode','auto','color','white');

%% What sequences went into experiment? Master table of barcodes, sequences, names in Collated.tsv:
x = readtable( 'Raw/Collated.txt' );
dx = table2cell( x );
all_barcodes= dx(:,4);
all_names = dx(:,1)
% OH THESE HAVE 1384 sequences, not all 233x.


%% Get names for barcodes that made it into Gun's ivstability file.
clear names_ivs;
for i = 1:length(barcodes_ivs);
    idx = find(strcmp(all_barcodes,barcodes_ivs{i}));
    names_ivs{i} = all_names{idx};
end

%%
mRNA_type = [];
for i = 1:length( names_ivs )
    name = upper(names_ivs{i});
    if ~isempty(strfind( name, 'MEV') )
        mRNA_type(i) = 1;
    elseif ~isempty(strfind( name, 'EGFP') )
        mRNA_type(i) = 5;
    elseif ~isempty(strfind( name, 'COV2') ) &&  ~isempty(strfind( name, 'F30') )
        mRNA_type(i) = 3;
    elseif ~isempty(strfind( name, 'F30') )
        mRNA_type(i) = 4;
    else
        assert( ~isempty(strfind(name,'NLUC')) );
        mRNA_type(i) = 2;
    end
end
mRNA_type_labels= {'MEV','Barna lab NLuc','COV2 Eterna','NLuc Eterna','eGFP'};

%% we should fit to exponentials and then order
clear tau_fit;
p0 = 6*60; % starter 
which_timepts = [1:6];
for i = 1:size(d_norm,1)
    tau_fit(i) = fminsearch( 'expfit', p0, [], times(which_timepts), d_norm(i,which_timepts) );
end

clf
plot( tau_fit/60, mean(d_norm(:,5:8)'),'.' )
xlabel( 'Time constant (hours)');
ylabel( 'Mean frac. intact, 3-6 hours' );
xlim([0 12]);
set(gcf, 'PaperPositionMode','auto','color','white');

%% Make heatmap of RNA's that made it through Gun's analysis
%   also we should separate by RNA groupings

set(figure(1),'position',[27     8   709   947]);
clf;
%[~,reorder] = sort( sum(d_norm(:,5:8)') ); % just persistence -- frac at end
%[~,reorder] = sort( tau_fit );
%[~,reorder] = sortrows([mRNA_type; tau_fit]')

set(groot,'defaultaxesticklabelinterpreter','none')
[~,reorder] = sortrows([mRNA_type; sum(d_norm(:,5:10)')]')
imagesc( log(d_norm(reorder,:))/log(2),[-5 5] );
set(gca,'ytick',1:size(d_norm,1),'yticklabel',names_ivs(reorder));
set(gca,'xtick',1:size(d_norm,2),'xticklabel',times/60);
xlabel('Degradation time (hours)')
set(gcf, 'PaperPositionMode','auto','color','white');
set(gca,'fontsize',6,'fontweight','bold');
set(gca,'position',[0.55 0.05 0.4 0.9]);

make_lines_horizontal(find( mRNA_type(reorder(1:end-1)) ~= mRNA_type(reorder(2:end)) ))

colormap(customcolormap_preset('red-white-blue'))
colorbar()


%% Let's take a closer look at some fits
set(figure(2),'pos',[200 200 800 500]);
idx = [];
colors = lines;
times_fine = [0:0.1:25]*60;
clf
designs = { 'hHBB_Nluc_hHBBx2', 'CoV2_P4_Nluc_hHBB','Curevac_2','Curevac_5', 'dG_max_3' };
designs = { 'Yellowstone', 'LinearDesign-1' };
%designs = { 'Curevac_2','Curevac_5','bpunpaired_min_2','9901251_2'};
%designs = { 'ModernaRef_egfp_lowstructure','ModernaRef_egfp_highstructure'};
clear titles
for q = 1:length(designs)
    m = find(strcmp( names_ivs, designs{q}));
    if isempty(m);    m = find(contains( names_ivs, designs{q}));   end
    if isempty(m); fprintf( [designs{q},'\n'] ); end;
    idx(q) = m(1);
    titles{q} = sprintf('%3.1f h: %s\n',tau_fit(idx(q))/60, names_ivs{idx(q)}  );
end
for q = 1:length(designs)
    plot( times/60, d_norm(idx(q),:),'o','markerfacecolor',colors(q,:),'markersize',12); hold on
end
for q = 1:length(designs)
    [~,pred] = expfit( tau_fit(idx(q)), times_fine );
    plot( times_fine/60,pred,'-','color',colors(q,:),'linew',2 );
end
h=legend(titles);set(h,'interpreter','none');
set(gca,'fontweight','bold','fontsize',9);
xlabel( 'Time (hours)' );
ylabel( 'Frac. intact (RT)');
set(gcf, 'PaperPositionMode','auto','color','white');


%% Do exponential fits.
set(figure(8),'pos',[62           1        1462         954]);
clf;
NBOOTSTRAP = 100;
semilogy( times/60, d_norm','o-' );
xlabel( 'Time (hours)' ); ylabel( 'Frac. full-length RT-able' );
set(gcf, 'PaperPositionMode','auto','color','white');
ylim([1e-3 1e1]);

%%
% we should fit to exponentials and then order
%clear tau_fit;
p0 = 6*60; % starter
which_timepts = [1:8];
for i = 1:size(d_norm,1)
    tau_fit(i) = fminsearch( 'expfit', p0, [], times(which_timepts), d_norm(i,which_timepts) );
    
    % bootstrap
    tau_boot = [];
    if NBOOTSTRAP > 0
        fprintf( 'Calculating %d bootstraps...\n',NBOOTSTRAP);
        for n = 1:NBOOTSTRAP
            boot_set = randi( length(which_timepts), length(which_timepts), 1 );
            tau_boot(n) = fminsearch( 'expfit', p0, [], times(which_timepts(boot_set)), d_norm(i,which_timepts(boot_set)) );
        end
        tau_fit_mean_boot(i) = mean(tau_boot);
        tau_fit_err_boot(i) = std(tau_boot);
    end
end
    
plot( tau_fit/60, mean(d_norm(:,5:8)'),'.' )
xlabel( 'Time constant (hours)');
ylabel( 'Mean frac. intact, 3-6 hours' );
xlim([0 12]);

%%
k_deg = 60./tau_fit;
k_deg_err = k_deg .* (tau_fit_err_boot./tau_fit);


%% Filter out 10 of 203 sequences that have anomalously high fraction intact at last timepoint,
%  1 sequence from COV2-Eterna study, and 1 sequence that is related to ENE Wilusz and caused confusion with barcoding
% Final table has 191 entries.
resid24h = d_norm(:,end) - exp( -times(end)./tau_fit' );
goodidx = find(resid24h < 0.05 & mRNA_type' ~= 3);

%% output data for paper table.
fprintf('\n')
%for i = 1:length(names_ivs)
filename = 'Rhiju_233x_stability_analysis_OUTPUT.txt';
fid = fopen( filename,'w' );
for i = goodidx'
    fprintf( '%d\t%s\t%s\t%f\t%f\n',i,barcodes_ivs{i},names_ivs{i},k_deg(i),k_deg_err(i));
    fprintf( fid, '%d\t%s\t%s\t%f\t%f\n',i,barcodes_ivs{i},names_ivs{i},k_deg(i),k_deg_err(i));
end
fclose(fid);
fprintf('Outputted %d rows to %s.\n',length(goodidx),filename)

