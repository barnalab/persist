function [loss,fit] = doubleexpfit(p,times,data);
% [loss,fit] = doubleexpfit(p,times,data);
%
% loss = rmse + 0.1 * abs(f_long)   
%
% abs(f_long) is L1 regularizer to keep fraction long low.

tau_short= p(1);
tau_long = p(2);
f_long   = p(3);
f_all   = p(4);

fit = f_all * ( (1-f_long) * exp( -times/tau_short ) + f_long * exp(-times/tau_long) );

loss = 0.0;
if ~exist( 'data','var') return; end;

rmse = sqrt( sum(( fit-data ).^2 ));
%rmse = sqrt( sum(( log(fit/data)/log(2) ).^2 ));
loss = rmse + 0.1 * abs(f_long);
