# Carrell2017
Matlab code from Carrell et al., 2017
Michael O'Connell 
08/02/2016

Contents:
dsat_2D_copy.m
dsat_jp2.m
dir: Other Files
results999999_dsat_jp2_mdo.mat
script_isres.m
sub_dsat2.sh

dsat_2D_v3.m 
The 2D version of the shuttling model that shows a shuttling phenotype. 

dsat_jp2.m
This file contains the nested-functions version of the shuttling (formerly “deconvolution 
saturated” or “dsat”) model, in which the
input parameters are assumed to be on a LOG scale. It also uses a Jacobian Pattern matrix, 
which tells the ODE solver the structure of the Jacobian. The published results are found
by running opening deconvolutionResults.mat and using xb(row,:) as input to 
decon_log_nest_JPattern(). 

Other Files: contains some files necessary to run the other .m files, such as gregsdata.mat
and hybrid_embryo.mat

results999999_dsat_jp2_mdo.mat


script_isres.m
To be used on the HPC, this script runs the ISRES optimization function on the shuttling
model (or any other model you want to feed it; use it as a template).
