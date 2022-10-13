% script_run_dsat
%
% This script starts running the function "run_dsat".


load('Mat/results999999_dsat_jp2_mdo.mat')

params = xb(1,1:20);
[error,penalty,dlwt,XT,beta,varnames,dl1x] = dsat_jp2(params);
plotComparison(dlwt.dlNuc,dlwt.dlCactNuc,XT);


params = 10.^xb(1,1:20);
Y = run_dsat(params,1,1);