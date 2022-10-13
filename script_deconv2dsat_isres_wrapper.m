% script_dsat_noshuttling
%
% This script will attempt to run evol opt code on the Toll sat model, but
% with the diffusivity of dl/Cact complex set to zero (actually,
% technically it's set to between 1e-12 and 1e-11).

yesplot = false;
close all
addpath('./Other Files')
addpath('C:\Users\gtreeves\Documents\Dropbox\Matlab\Dorsal\MDO\deconvolution')
addpath('C:\Users\gtreeves\Documents\Dropbox\Matlab\Dorsal\MDO\deconvolution\Other Files')

%
% Load parameters
%
load('Mat/2017-02-10_EP_deconv')
load('deconvolutionResults.mat')
[~,isort] = sort(E(1:i-1));
N = i - 1;
nparams = round(0.05*N);
XB = xb;

for ii = 15:nparams
	params = XB(ii,1:end-1);
	script_deconv2dsat_isres
end