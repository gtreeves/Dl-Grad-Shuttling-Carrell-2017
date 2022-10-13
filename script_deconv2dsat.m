	% script_deconv2dsat
%
% This script takes a deconvolution model paramter set and converts it to a
% dsat set.  There is some futz'ing that has to be done with adding the
% four extra parameters.  To make it happen, we hold the original 15
% parameters constant, and start with omega and rho large (like, 1000?).
% The parameters are also subject to the constraint that Vc =
% Atoll*epsilon*eta/rho.  This constraint cannot be perfectly satisfied, as
% both Vc and Atoll change in time in an uncorrelated way, but the
% parameter values can probably be chosen to *almost* satisfy the
% constraint.
%
% nu is completely unconstrained, so we will start with nu = 1.  We will
% also set epsilon = eta.
%
% Once we start with this parameter set, we will allow omega and rho to
% become meso.  We will see what happens then.

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

%
% Load constraints on eta, epsilon. 
%
% Vc = Acs*epsilon*eta/rho
%
load('Mat/nuclearsize_vars')
Vc2 = cell2mat(struct2cell(Vc));
Acs2 = cell2mat(struct2cell(Acs));
rho0 = 50; omega0 = rho0;
eps0 = sqrt(Vc2(end)/Acs2(end)*rho0);
eta0 = eps0;
nu0 = 1;
extras = log10([nu0 omega0 eta0 rho0 eps0]);
n_extras = length(extras);

xbd = xb(isort,:);
xbd = xbd(1:nparams,:);

nSteps = 20;
Extra_params = cell(nparams,1);
E_set = zeros(nparams,nSteps);
P2_set = zeros(nparams,nSteps);
for ii = 1:nparams
	params = xbd(ii,1:end-1);
	
	
% 	figure('pos',1.0e+03 * [0.1546    0.0474    1.1792    0.7264])
	
	%
	% Deconv model, then dsat model
	%
	[e_deconv,p_deconv,dlwtd,XTd,betad,varnamesd,dl1xd] = decon_log_nest_JPattern(params);
	Params = [params extras];
	[e_dsat,p_dsat,dlwt,XT,beta,varnames,dl1x] = dsat_jp2(Params);
	
	%
	% Plotting, if asked for.
	%
	if yesplot
		error1 = plotComparison(dlwtd.dlNuc,dlwtd.dlCactNuc,XTd);
		error1 = plotComparison(dlwt.dlNuc,dlwt.dlCactNuc,XT);
		
		t = XTd.T14;
		xd = XTd.X14(1,:)';
		
		nt = length(t);
		i14 = round(nt*0.99);
		y2xd = dlwtd.dlNuc.NC14(:,i14);
		y1xd = dl1xd.dlNuc.NC14(:,i14);
		
		t = XT.T14;
		x = XT.X14(1,:)';
		
		nt = length(t);
		i14 = round(nt*0.99);
		y2x = dlwt.dlNuc.NC14(:,i14);
		y1x = dl1x.dlNuc.NC14(:,i14);
		
		figure('pos',[154.6000   47.4000  483.2000  726.4000])
		% 	figure('pos',[154.6000   47.4000  1231.2  726.4000])
		subplot(2,1,1)
		plot(xd,y2xd,x,y2x)
		title('Toll sat, no shuttling')
		hold on
		plot(xd,y1xd,x,y1x)
		legend('deconv 2x','dsat 2x','deconv 1x','dsat 1x')
		
		subplot(2,1,2)
		plot(xd,mDU(y2xd),x,mDU(y2x))
		title('Normalized')
		hold on
		plot(xd,mDU(y2xd),x,mDU(y1x))
		legend('deconv 2x','dsat 2x','deconv 1x','dsat 1x')
	end
	
	%
	% Now that we have the starting params, we need to alter them to go
	% towards having Toll saturate more, so that dosage affects it, without
	% changing the error.  We also want the second element of the penalty
	% function to be ever-decreasing.
	%
	% rho0 = 100000; omega0 = rho0;
	% eps0 = sqrt(Vc2(end)/Acs2(end)*rho0);
	% eta0 = eps0;
	% nu0 = 1;
	% extras = log10([nu0 omega0 eta0 rho0 eps0]);
	%
	extras0 = extras;
	delt = -0.1;
	delt1 = 1e-4;
	rho_ratio = 1 + delt;
	param_ratio = 1 + delt1;
	
	Extras = zeros(nSteps,n_extras);
	E = zeros(nSteps,1);
	P2 = zeros(nSteps,1);
	Extras(1,:) = extras; E(1) = e_dsat; P2(1) = p_dsat(2);
	for i = 2:nSteps
		
		%
		% First, get the gradient vectors for the error (g) and the dosage
		% penalty (h).
		%
		g = zeros(1,n_extras);
		h = g;
		for k = 1:n_extras
			extras1 = extras;
			extras1(k) = extras1(k) + log10(param_ratio);
			
			Params = [params extras1];
			[e_dsat1,p_dsat1] = dsat_jp2(Params);
			
			g(k) = (e_dsat1 - e_dsat)/(10^extras(k)*delt1);
			h(k) = (p_dsat1(2) - p_dsat(2))/(10^extras(k)*delt1);
		end
		
		%
		% Now, systematically decrease rho and omega.
		%
		extras(2) = extras(2) + log10(rho_ratio);
		extras(4) = extras(4) + log10(rho_ratio);
		
		%
		% Now find how to change the other parameters while keeping the
		% error the same and also decreasing the dosage penalty.  We will
		% use the steepest descent in penalty with the constraint of
		% holding error constant.
		%
		de = -(g(2)*10^extras(2) + g(4)*10^extras(4))*delt;
		H = g([1 3 5])*h([1 3 5])';
		a = -de/H;
		
		d_param = -a*h;
		if any(d_param < -1)
			b = -1./d_param;
			b = min(b(b > 0));
			
			d_param = d_param*b/2;
		end
		
		d_logparam = log10(1 + d_param([1 3 5])./10.^extras([1 3 5]));
		extras([1 3 5]) = extras([1 3 5]) + d_logparam;
		
		Params = [params extras];
		[e_dsat1,p_dsat] = dsat_jp2(Params);
		E(i) = e_dsat1;
		P2(i) = p_dsat(2);
		
		Extras(i,:) = extras;
		
	end
	Extra_params{ii} = Extras;
	E_set(ii,:) = E';
	P2_set(ii,:) = P2';
	
	%
	% Plot a given parameter set?
	%
	yesplot1 = yesplot;
	yesplot = false;
	if yesplot
		Extras = Extra_params{5};
		Params = [params Extras(7,:)];
		[e_dsat,p_dsat,dlwt,XT,beta,varnames,dl1x] = dsat_jp2(Params);
		error1 = plotComparison(dlwt.dlNuc,dlwt.dlCactNuc,XT);
		
	end
	yesplot = yesplot1;

	disp(ii)
end


% NC = {'10','11','12','13','14'};
% figure
% for i = 1:length(NC)
% nc1 = NC{i};
% 
% nc = ['NC',nc1];
% xc = ['X',nc1];
% x = XT.(xc)(1,:);
% 
% Y = dlwt.tollAct;
% yy = Y.(nc);
% 
% X = dlwt.tollDC;
% xx = X.(nc);
% 
% W = dlwt.cactCyt;
% wc = W.(nc);
% 
% nu = 10^(Params(16));
% omega = 10^(Params(17));
% eta = 10^(Params(18));
% 
% plot(x,eta*yy.*wc-(nu+omega)*xx)
% hold on
% end
