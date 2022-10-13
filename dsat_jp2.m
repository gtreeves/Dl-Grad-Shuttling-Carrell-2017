function [ varargout ] = dsat_jp2(Params,varargin)
% DSAT_JP2 Deconvolution model w/ saturdated Toll receptors using JPattern
% to specify the structure of the Jacobian, version 2. 
%
%   [ERROR,PENALTY] = DSAT_JP2(PARAMS) takes a vector of
%   parameters (PARAMS), where size(PARAMS)=[1 20], and calculates the RSS
%   error (ERROR) against those found in GREGSDATA; the PENALTY
%   indicates if the solution is feasible/infeasible, respectively.
%
%   [ERROR,PENALTY,PROTEIN2X,TIME,BETA,NAMES,PROTEIN1X] = DSAT_JP2(PARAMS) 
%   creates nested structures PROTEIN2X and PROTEIN1X, in which each species and
%   nuclear cycle is accessible as PROTEIN.SPECIES.NC. The list of species
%   is 'dlNuc','dlCyt','dlCactNuc','dlCactCyt','cactNuc','cactCyt',
%   'tollAct', and 'tollDC'. Each species has NC fields 'NC10','NC10m',
%   'NC11','NC11m','NC12','NC12m','NC13','NC13m',and 'NC14'.
%   2X and 1X refer to 2x dl (wt) and 1x dl genotypes. TIME is a similar
%   structure that refers to the accompanying X and T values (spacial
%   coordinates and time vectors) for each species & nuclear cycle. BETA
%   is the value that scales the simulation output to the microscopy data
%   from GREGSDATA. NAMES is a cell variable that contains the names of the
%   parameters in the model.
%
% GTR change on 2/16/2017: changed the penalty function so that the dosage
% penalty gets stored in penalty(2) instead of adding to penalty(1).

%% misc. initial procedures

if length(Params) > 20
	dl0 = Params(21);
	Params(21:end) = [];
else
	dl0 = 1;
end

addpath('C:\Users\gtreeves\Documents\Dropbox\Matlab\Dorsal\MDO\shuttling/Other Files')
opts = optimset('Display','off');
args = struct('dl', dl0, 'cact',1,'genefun',[],'plot',false);
Beta = 1;
if nargin > 1
    if mod(length(varargin),2)==0
        for i = 1:2:length(varargin)
            args.(varargin{i})=varargin{i+1};
        end
    else
        error('Varargin requires name-value pairs. Some options not specified')
    end
end
options = [];
varargout = cell(1,nargout);
E = zeros(size(Params,1),1); penalty = zeros(size(Params,1),3);
% Greg's data
[C1,dC1] = gregsData('nuclear');
[C2,dC2] = gregsData('total');
% C = C1; dC = dC1;
C = [C1;C2]; dC = [dC1;dC2]; 
% clear C1 C2 dC1 dC2


ncs = {'NC10','NC10m','NC11','NC11m','NC12','NC12m','NC13','NC13m','NC14'};
if nargout <= 2 && ~args.plot
    names = {'dlNuc','dlCyt','dlCactNuc'};
    times = false;
    
else
    names = {'dlNuc','dlCyt','dlCactNuc','dlCactCyt','cactNuc','cactCyt',...
        'tollAct','tollDC'};
    ts = {'T10','T10m','T11','T11m','T12','T12m','T13','T13m','T14'};
    xs = {'X10','X10m','X11','X11m','X12','X12m','X13','X13m','X14'};
    nt=1;
    times = true;
end
% options = odeset('RelTol',1e-1,'AbsTol',1e-1);

load('Mat/hybrid_embryo.mat','data')
nc = ncs{1};
% Constants
% points in x (# of nuclei)
M = struct('NC10',13,'NC10m',13,'NC11',19,'NC11m',19,'NC12',26,'NC12m',...
    26,'NC13',36,'NC13m',36,'NC14',51);
% Length (L2 - L1) in x
L1 = 0; L2 = 1;
% time spans
tspan = struct('NC10',data.t(1:16),'NC10m',data.t(16:33),'NC11',data.t(1:16)+data.t(33),...
    'NC11m',data.t(17:33)+data.t(33),'NC12',data.t(34:59)+data.t(33),...
    'NC12m',data.t(60:76)+data.t(33),'NC13',data.t(77:126)+data.t(33),...
    'NC13m',data.t(127:148)+data.t(33),'NC14',linspace(35.6,90,184)'+data.t(33));
clear data

%initialize protein
sub_pro = struct('NC10',zeros(13,16),'NC10m',zeros(13,18),'NC11',...
    zeros(19,16),'NC11m',zeros(19,17),'NC12',zeros(26),'NC12m',...
    zeros(26,17),'NC13',zeros(36,50),'NC13m',zeros(36,22),'NC14',...
    zeros(51,184));
protein = struct('dlNuc',sub_pro,'dlCyt',sub_pro,'dlCactNuc',sub_pro,...
    'dlCactCyt',sub_pro,'cactNuc',sub_pro,'cactCyt',sub_pro,...
    'tollAct',sub_pro,'tollDC',sub_pro,'nuclearDorsal',sub_pro);
protein2x = protein;
clear sub_pro

%initialize 
time = struct; 

% other nested variables
P=[];An = struct('NC10',[],'NC10m',[],'NC11',...
    [],'NC11m',[],'NC12',[],'NC12m', [],'NC13',[],'NC13m',[],'NC14', []);
Am=An;
Acs=An;
Vn=An;
Vc=An;
for i = 1:9
    if mod(i,2)==0
[An.(ncs{i}),Am.(ncs{i}),Acs.(ncs{i}),Vn.(ncs{i}),Vc.(ncs{i})]=nuclearSize(1,'static',M.(ncs{i}),'mitosis');
    else
[An.(ncs{i}),Am.(ncs{i}),Acs.(ncs{i}),Vn.(ncs{i}),Vc.(ncs{i})]=nuclearSize(1,'static',M.(ncs{i}),'interphase');
    end
end
%% MAIN
if ~isempty(Params)
for rows = 1:size(Params,1)
    try
%     if rows > 1
%         disp([E(rows-1) penalty(rows-1,:)])
%     end
    Params_ = 10.^(Params(rows,:));
    Params_ = num2cell(Params_);
    [lambdaU,lambdaW,lambdaV,sigmaU,sigmaW,sigmaV,muU,muW,muV,...
        gamma,psi,alpha,beta,phi,beta0,nu,omega,eta,rho,epsilon]=Params_{:};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2x Dorsal %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x = linspace(0,1,M.NC10)';
%     [An,Am,Acs,Vn,Vc]=nuclearSize(1,'static',M.NC10,'interphase');
    y0 = [zeros(M.NC10,1);
        zeros(M.NC10,1);
        args.dl*ones(M.NC10,1);
        args.dl*ones(M.NC10,1);
        args.cact*(1/alpha)*ones(M.NC10,1);
        args.cact*(1/alpha)*ones(M.NC10,1);
        (beta/rho).*Vc.NC10*rho/(Acs.NC10*epsilon*eta)*exp(-0.5*(x/phi).^2);
        zeros(M.NC10,1)];
    
    
	for i = 1:9
		if mod(i,2)==0
			nc = ncs{i};
			mitosis;
		else
			nc = ncs{i};
			interphase;
		end
	end
    
    for j=1:length(ncs)
        protein.('nuclearDorsal').(ncs{j})=protein.dlNuc.(ncs{j})+...
            protein.dlCactNuc.(ncs{j});
        protein.('totalDorsal').(ncs{j})...
            = (Vn.(ncs{j})*(protein.dlNuc.(ncs{j})+...
            protein.dlCactNuc.(ncs{j}))+Vc.(ncs{j})*(protein.dlCyt.(ncs{j})+...
            protein.dlCactCyt.(ncs{j})))/(Vc.(ncs{j})+Vn.(ncs{j}));
    end
    protein2x=protein;
    try
        % Simulation data as a column vector
        Dnuc = [protein.dlNuc.NC11(:);
            protein.dlNuc.NC12(:);
            protein.dlNuc.NC13(:);
            protein.dlNuc.NC14(:)];
        Dcyt = [protein.dlCactNuc.NC11(:);
            protein.dlCactNuc.NC12(:);
            protein.dlCactNuc.NC13(:);
            protein.dlCactNuc.NC14(:)];
        D = [protein.nuclearDorsal.NC11(:);
            protein.nuclearDorsal.NC12(:);
            protein.nuclearDorsal.NC13(:);
            protein.nuclearDorsal.NC14(:);
            protein.totalDorsal.NC11(:);
            protein.totalDorsal.NC12(:);
            protein.totalDorsal.NC13(:);
            protein.totalDorsal.NC14(:)];

        
        Beta = mean(D.*C)/mean(D.^2);
        E(rows) = sqrt(sum(((C-Beta*D)./dC).^2)/length(D));
        
        % if dlCactNuc is closer to gregsData('nuclear') than dlNuc, penalize.
        Beta1 = mean(Dnuc.*C1)/mean(Dnuc.^2);
        Beta2 = mean(Dcyt.*C1)/mean(Dcyt.^2);
        penalty(rows,1) = penalty(rows,1) + ...
            (sqrt(sum(((C1-Beta1*Dnuc)./dC1).^2)/length(Dnuc)) > ...
            sqrt(sum(((C1-Beta2*Dcyt)./dC1).^2)/length(Dcyt)));
        
		%
		% Fit to a repressive Hill function, then get the half-max of the
		% fit. That is so we can compare this half-max to the analog for
		% the dorsal 1x situation.  The expansion in the 1x must be greater
		% than 1.5 nuclei (ie, x = 0.03) for there to be no penalty.
		%
        dg=protein.nuclearDorsal.NC14(:,end);
        [AA,BB,CC,DD]=fitstep(x,dg);
        dgfit = AA.*10.^(BB)./(10.^(BB)+x.^(2*CC))+DD;
        hmax = min(dgfit)+0.5*(max(dgfit)-min(dgfit));
        sft = @(x) AA.*10.^(BB)./(10.^(BB)+x.^(2*CC))+DD-hmax;
        zed2=fzero(sft,0.5,opts);
    catch ME
        
        if strcmp(ME.message,'Matrix dimensions must agree.') || strcmp(ME.message,'Exiting due to infeasibility: 1 lower bound exceeds the corresponding upper bound.')
        E(rows) = 1e6+randn;
        penalty(rows,:)=1e6+randn;
        else
            error(ME.message)
        end
        continue
    end
%     if args.plot && rows == 1 && times
%         protNames = fieldnames(protein);
%         for k = 1:length(protNames)
% %             figure; hold on;
%             for j = 1:9
%                 surf(time.(ts{j})',time.(xs{j})',Beta*protein.(protNames{k}).(ncs{j}))
%                 surf(time.(ts{j})',-time.(xs{j})',Beta*protein.(protNames{k}).(ncs{j}))
%             end
%             shading flat
%             view([45 45])
%             title(protNames{k})
%         end
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1x Dorsal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = linspace(0,1,M.NC10)';
%     [An,Am,Acs,Vn,Vc]=nuclearSize(1,'static',M.NC10,'interphase');
    y0 = [zeros(M.NC10,1);
        zeros(M.NC10,1);
        args.dl*0.5*ones(M.NC10,1);
        args.dl*0.5*ones(M.NC10,1);
        args.cact*(1/alpha)*ones(M.NC10,1);
        args.cact*(1/alpha)*ones(M.NC10,1);
        (beta/rho).*Vc.NC10*rho/(Acs.NC10*epsilon*eta)*exp(-0.5*(x/phi).^2);
        zeros(M.NC10,1)];
    
    nt=1;
    for i = 1:9                                                         
        if mod(i,2)==0                                                  
            nc = ncs{i};                                                
            mitosis;                                                    
        else                                                            
            nc = ncs{i};                                                
            interphase;                                                 
        end                                                             
    end                                                                 
        
    for j=1:length(ncs)
        protein.('nuclearDorsal').(ncs{j})=protein.dlNuc.(ncs{j})+...
            protein.dlCactNuc.(ncs{j});
        protein.('totalDorsal').(ncs{j})...
            = (Vn.(ncs{j})*(protein.dlNuc.(ncs{j})+...
            protein.dlCactNuc.(ncs{j}))+Vc.(ncs{j})*(protein.dlCyt.(ncs{j})+...
            protein.dlCactCyt.(ncs{j})))/(Vc.(ncs{j})+Vn.(ncs{j}));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Penalty
 
%     vtot = protein2x.nuclearDorsal.NC14(1,:); vtot = vtot-min(vtot); vtot = vtot/max(vtot);
%     vdl = protein2x.dlNuc.NC14(1,:); vdl = vdl-min(vdl); vdl = vdl/max(vdl);
%     
%     btot = protein2x.nuclearDorsal.NC14(end,:); btot = btot-min(btot); btot = btot/max(btot);
%     bdl = protein2x.dlCactNuc.NC14(end,:); bdl = bdl-min(bdl); bdl = bdl/max(bdl);
% 
%     penalty(rows,1) = sqrt((sum((vtot-vdl).^2)+sum((btot-bdl).^2)))-...
%         sqrt((sum((vtot-bdl).^2)+sum((btot-vdl).^2)));
    dg=protein.nuclearDorsal.NC14(:,end);
    try
		
	%
	% Fit to a repressive Hill function, then get the half-max of the fit.
	% Compare this half-max to the analog for the dorsal 2x situation.  The
	% expansion in the 1x must be greater than 1.5 nuclei (ie, x = 0.03)
	% for there to be no penalty.
	%
    [AA,BB,CC,DD]=fitstep(x,dg);
    dgfit = AA.*10.^(BB)./(10.^(BB)+x.^(2*CC))+DD;
    hmax = min(dgfit)+0.5*(max(dgfit)-min(dgfit));
    sft = @(x) AA.*10.^(BB)./(10.^(BB)+x.^(2*CC))+DD-hmax;
    zed1=fzero(sft,0.5,opts);
    z = zed2-zed1;
    catch ME
        z = NaN;
        if ~strcmp(ME.message,'Exiting due to infeasibility: 1 lower bound exceeds the corresponding upper bound.')
            error(ME.message)
        end
    end
    if isnan(z)
        z = 1;
    end
%     if E(rows) < 3 % wait until we have a decent result before calculating this.
        penalty(rows,2) = penalty(rows,2)+ z+0.03; % adding a small number to force the expansion
%     end
    catch ME2
        disp(ME2.message)
        E(rows) = NaN;
        penalty(rows,1)=5;
%         penalty(rows,1)=5;
    end
%     penalty(rows,2) = 0;
    
    if protein2x.dlNuc.NC14(end,end) > protein2x.dlCactNuc.NC14(end,end)
        penalty(rows,3) = 1;
    end
end
end

%% Gene expression
% if ~isempty(args.genefun)
%     dlgrad = protein2x.dlNuc.NC14(:,end);
%     load('hybrid_embryo.mat','data')
% 
%   genename = {'sna', 'sog', 'dpp', 'vnd'};
% % s0 = [0 0.3 1 0.27]; % from Greg's code
% % initialize values for loop
% geneavg = zeros(151,4);
% 
% for i = 1:4
%     matname = [genename{i},'avg.mat'];
%     gg=load(matname,'t');
%     if i == 2
%         geneavg(1:121,i) = gg.t(31:end);
%     else
%         geneavg(:,i) = gg.t;
%     end
%     
% end
% geneavg(1:25,1)=round(geneavg(1:25,1));
% geneavg(133:151,3)=round(geneavg(133:151,3));
% geneavg = interp1(linspace(0,1,151),geneavg,linspace(0,1,51));
% %     gp =[ thetaDlSna,thetaDlSog,thetaSnaSog,thetaDlDpp,...
% %         thetaSnaVnd,thetaDlVnd,tauSna,tauSog,tauDpp,tauVnd,...
% %         nSna,nSog,nDpp,nVnd,noise];
%     gp0 = log10([dlgrad(11) dlgrad(24) 0.5 dlgrad(31) 0.5 dlgrad(18) 0.1,1,1,1,100,100,100,100,0.2]);
%     fh = @(x) genef_log_nest_avg(gp0,protein2x.dlNuc);
%     [gp,fval]=fmincon(fh,gp0,[],[],[],[],-3*ones(1,15),ones(1,15));
%     [errg,peng,proteinsg]=genef_log_nest_avg( gp,protein2x.dlNuc );
% end

%% Output

pnames = {'lambdaU','lambdaW','lambdaV','sigmaU','sigmaW','sigmaV',...
    'muU','muW','muV','gamma','psi','alpha','beta','phi','beta0','nu',...
    'omega','eta','rho','epsilon'};
argsout = {E,penalty,protein2x,time,Beta,pnames,protein};

for i = 1:nargout
    varargout{i} = argsout{i};
end

%% nested functions
    function interphase
        %update variables
        x = linspace(L1,L2,M.(nc));
        e = ones(M.(nc),1);P = spdiags([e -2*e e],[-1 0 1],M.(nc),M.(nc));
        j1 = spdiags(ones(M.(nc),1),0,M.(nc),M.(nc));
        options = odeset('RelTol',1e-1,'AbsTol',1e-2,'JPattern', [P j1 j1 j1 j1 j1 j1 j1;
            j1 P j1 j1 j1 j1 j1 j1;
            j1 j1 P j1 j1 j1 j1 j1;
            j1 j1 j1 P j1 j1 j1 j1;
            j1 j1 j1 j1 P j1 j1 j1;
            j1 j1 j1 j1 j1 P j1 j1;
            j1 j1 j1 j1 j1 j1 P j1;
            j1 j1 j1 j1 j1 j1 j1 P]);
        
        P(1,2) = 2;P(M.(nc),M.(nc)-1) = 2;
        
%         [An.(ncs{i}),Am.(ncs{i}),Acs.(ncs{i}),Vn.(ncs{i}),Vc.(ncs{i})]=nuclearSize(1,'static',M.(nc),'interphase');
        
        % interpolating (nuc = cyt)
        y0_ = reshape(y0,length(y0)/8,8);
        y0 = zeros(M.(nc)*8,1);
        
        for n = 1:2:6
            y0((n-1)*M.(nc)+1:n*M.(nc),1)=interp1(linspace(0,1,...
                size(y0_,1)),y0_(:,n+1),linspace(0,1,M.(nc)));
            y0((n)*M.(nc)+1:(n+1)*M.(nc),1)=interp1(linspace(0,1,...
                size(y0_,1)),y0_(:,n+1),linspace(0,1,M.(nc)));
        end
        for n = 7:8
            y0((n-1)*M.(nc)+1:n*M.(nc),1)=interp1(linspace(0,1,...
                size(y0_,1)),y0_(:,n),linspace(0,1,M.(nc)));
        end
        
        if times
            [t,y] = ode15s(@int_fun,tspan.(nc),y0,options);
            [time.(xs{nt}),time.(ts{nt})] = meshgrid(x,t);
            nt=nt+1;
        else
            [~,y] = ode15s(@int_fun,tspan.(nc),y0,options);
        end
        
        y = y';
        
        for i_ = 1:length(names)
            protein.(names{i_}).(nc)=y((i_-1)*M.(nc)+1:(i_)*M.(nc),:);
        end
        
        % update y0 for mitosis (mixing)
        y0 = y(:,end);
        for i_ =[6 5 4 3 2 1]
            if mod(i_,2)==0
                y0((i_-1)*M.(nc)+1:(i_)*M.(nc)) = ...
                    (Vn.(nc)*y0((i_-2)*M.(nc)+1:(i_-1)*M.(nc)) +...
                    Vc.(nc)*y0((i_-1)*M.(nc)+1:(i_)*M.(nc)))/(Vn.(nc)+Vc.(nc));
            else
                y0((i_-1)*M.(nc)+1:(i_)*M.(nc))=zeros(M.(nc),1);
            end
        end
        
    end
    function F = int_fun(~,y)
        
        un = y(1:M.(nc));
        uc = y(M.(nc)+1:2*M.(nc));
        wn = y(2*M.(nc)+1:3*M.(nc));
        wc = y(3*M.(nc)+1:4*M.(nc));
        vn = y(4*M.(nc)+1:5*M.(nc));
        vc = y(5*M.(nc)+1:6*M.(nc));
        yy = y(6*M.(nc)+1:7*M.(nc));
        xx = y(7*M.(nc)+1:8*M.(nc));
        
        f1 = (An.(nc)*(sigmaU*uc-muU*un)+Vn.(nc)*(beta0*wn-gamma*un.*vn))/(Vn.(nc));
        f2 = (lambdaU*Am.(nc)*P*uc+omega*epsilon*Acs.(nc)*xx+Vc.(nc)*(beta0*wc-gamma*uc.*vc)...
            -An.(nc)*(sigmaU*uc-muU*un))/(Vc.(nc));
        f3 = (An.(nc)*(sigmaW*wc-muW*wn)-Vn.(nc)*(beta0*wn-gamma*un.*vn))/(Vn.(nc));
        f4 = (lambdaW*Am.(nc)*P*wc-Acs.(nc)*epsilon*(eta*wc.*yy-nu*xx)...
            -Vc.(nc)*(beta0*wc-gamma*uc.*vc)-An.(nc)*(sigmaW*wc-muW*wn))/(Vc.(nc));
        f5 = (An.(nc)*(sigmaV*vc-muV*vn)+Vn.(nc)*psi*(beta0*wn-gamma*un.*vn))/(Vn.(nc));
        f6 = (lambdaV*Am.(nc)*P*vc+psi*(omega*epsilon*Acs.(nc)*xx+Vc.(nc)*(beta0*wc-gamma*uc.*vc))+1-alpha*Vc.(nc)*vc...
            -An.(nc)*(sigmaV*vc-muV*vn))/(Vc.(nc));
        f7 = (nu+omega)*xx-eta*yy.*wc+beta*exp(-0.5*(x'/phi).^2)-rho*yy;
        f8 = eta*yy.*wc-(nu+omega)*xx;
        
        F = [f1;f2;f3;f4;f5;f6;f7;f8];
    end

    function mitosis
        % update variables
%         [An.(ncs{i}),Am.(ncs{i}),Acs.(ncs{i}),Vn.(ncs{i}),Vc.(ncs{i})]=nuclearSize(1,'static',M.(nc),'mitosis');
        %         options = odeset('Jacobian',@mit_jac);
        if times
            [t,y] = ode15s(@mit_fun,tspan.(nc),y0,options);
            [time.(xs{nt}),time.(ts{nt})] = meshgrid(x,t);
            nt=nt+1;
        else
            [~,y] = ode15s(@mit_fun,tspan.(nc),y0,options);
        end
        y = y';
        
        for i_ = 1:length(names)
            protein.(names{i_}).(nc)=y((i_-1)*M.(nc)+1:(i_)*M.(nc),:);
        end
        
        % update y0 for interphase (interpolation)
    end

    function F = mit_fun(~,y)
        
        uc = y(M.(nc)+1:2*M.(nc));
        wc = y(3*M.(nc)+1:4*M.(nc));
        vc = y(5*M.(nc)+1:6*M.(nc));
        yy = y(6*M.(nc)+1:7*M.(nc));
        xx = y(7*M.(nc)+1:8*M.(nc));
        
        
        f1 = zeros(M.(nc),1);
        f2 = (lambdaU*Am.(nc)*P*uc+omega*epsilon*Acs.(nc)*xx-gamma*Vc.(nc)*uc.*vc+beta0*Vc.(nc)*wc)/(Vc.(nc));
        f3 = zeros(M.(nc),1);
        f4 = (lambdaW*Am.(nc)*P*wc-eta*epsilon*Acs.(nc)*wc.*yy+nu*epsilon*Acs.(nc)*xx+gamma*Vc.(nc)*uc.*vc-beta0*Vc.(nc)*wc)/(Vc.(nc));
        f5 = zeros(M.(nc),1);
        f6 = (lambdaV*Am.(nc)*P*vc+psi*(omega*epsilon*Acs.(nc)*xx-gamma*Vc.(nc)*uc.*vc+beta0*Vc.(nc)*wc)+1-...
            alpha*Vc.(nc)*vc)/(Vc.(nc));
        f7 = (nu+omega)*xx-eta*yy.*wc+beta*exp(-0.5*(x'/phi).^2)-rho*yy;
        f8 = eta*yy.*wc-(nu+omega)*xx;
        
        F = [f1;f2;f3;f4;f5;f6;f7;f8];
        
    end
end