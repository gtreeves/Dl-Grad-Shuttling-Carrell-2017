function [ varargout ] = dsat_2D_copy
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
%   structure that refers to the accompanying Y and T values (spacial
%   coordinates and time vectors) for each species & nuclear cycle. BETA
%   is the value that scales the simulation output to the microscopy data
%   from GREGSDATA. NAMES is a cell variable that contains the names of the
%   parameters in the model.

Params = [-2.73893196930974,2.94903962777468,-0.603866134315929,-0.132948063685535,-0.834026108393399,0.883069371497378,-0.0866451744284050,-2.76139305374563,0.420647032100406,-1.16086197373076,2.07514898377416,0.0129872996158352,2.98662493032353,-0.875586131020344,0.165457783682081,-2.89345687083553,-1.16922779768909,2.57686028464904,-1.78047838373452,-2.54449321699653,2.98662493032353,-0.675586131020344];
sna_threshold = 1.11;
% Params(1:3) = 0;
%% misc. initial procedures
opts = optimset('Display','off');
args = struct('dl', 1, 'cact',1,'genefun',[],'plot',true);
Beta = 1;
% if nargin > 2
%     if mod(length(varargin),2)==0
%         for i = 1:2:length(varargin)
%             args.(varargin{i})=varargin{i+1};
%         end
%     else
%         error('Varargin requires name-value pairs. Some options not specified')
%     end
% end
options = [];
varargout = cell(1,nargout);
E = zeros(size(Params,1),1); penalty = E;



ncs = {'NC10','NC10m','NC11','NC11m','NC12','NC12m','NC13','NC13m','NC14'};
if nargout <= 2 && ~args.plot
    names = {'dlNuc','dlCyt','dlCactNuc'};
    times = false;
    
else
    names = {'dlNuc','dlCyt','dlCactNuc','dlCactCyt','cactNuc','cactCyt',...
        'tollAct','tollDC'};
    ts = {'T10','T10m','T11','T11m','T12','T12m','T13','T13m','T14'};
    xs = {'X10','X10m','X11','X11m','X12','X12m','X13','X13m','X14'};
    ys = {'Y10','Y10m','Y11','Y11m','Y12','Y12m','Y13','Y13m','Y14'};
    nt=1;
    times = true;
end
% options = odeset('RelTol',1e-1,'AbsTol',1e-1);

load('hybrid_embryo.mat','data')
nc = ncs{1};
% Constants
% points in x (# of nuclei)
M = struct('NC10',13,'NC10m',13,'NC11',19,'NC11m',19,'NC12',26,'NC12m',...
    26,'NC13',36,'NC13m',36,'NC14',51);
% N = struct('NC10',18,'NC10m',18,'NC11',26,'NC11m',26,'NC12',36,'NC12m',...
%     36,'NC13',50,'NC13m',50,'NC14',70);
% M = struct('NC10',5,'NC10m',5,'NC11',8,...
%     'NC11m',8,'NC12',12,'NC12m',12,'NC13',17,'NC13m',17,'NC14',25);
% N=M;
% N = struct('NC10',floor(18/2),'NC10m',floor(18/2),'NC11',floor(26/2),...
%     'NC11m',floor(26/2),'NC12',floor(36/2),'NC12m',floor(36/2),'NC13',...
%     floor(50/2),'NC13m',floor(50/2),'NC14',floor(70/2));
% M = struct('NC10',20,'NC10m',20,'NC11',20,'NC11m',20,'NC12',20,'NC12m',...
%     20,'NC13',20,'NC13m',20,'NC14',20);
% N=struct('NC10',40,'NC10m',40,'NC11',40,'NC11m',40,'NC12',40,'NC12m',...
%     40,'NC13',40,'NC13m',40,'NC14',40);
N=M;
% Length (L2 - L1) in x
% Greg's data
[C,dC] = gregsData('nuclear','M',[M.NC11 M.NC12 M.NC13 M.NC14]);
% [C2,dC2] = gregsData('total');
% C = C1; dC = dC1;
% C = [C1;C2]; dC = [dC1;dC2]; 
% clear C1 C2 dC1 dC2

% L1 = 0; L2 = 1; 
% time spans
tspan = struct('NC10',data.t(1:16),'NC10m',data.t(16:33),'NC11',data.t(1:16)+data.t(33),...
    'NC11m',data.t(17:33)+data.t(33),'NC12',data.t(34:59)+data.t(33),...
    'NC12m',data.t(60:76)+data.t(33),'NC13',data.t(77:126)+data.t(33),...
    'NC13m',data.t(127:148)+data.t(33),'NC14',linspace(35.6,90,184)'+data.t(33));
clear data

%initialize protein
sub_pro = struct('NC10',zeros(5,5,16),'NC10m',zeros(5,5,18),'NC11',...
    zeros(8,8,16),'NC11m',zeros(8,8,17),'NC12',zeros(12,12,26),'NC12m',...
    zeros(12,12,17),'NC13',zeros(17,17,50),'NC13m',zeros(17,17,22),'NC14',...
    zeros(25,25,184));
% sub_pro = struct('NC10',[],'NC10m',[],'NC11',...
%     [],'NC11m',[],'NC12',[],'NC12m',[],'NC13',[],'NC13m',[],'NC14',[]);
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
[An.(ncs{i}),Am.(ncs{i}),Acs.(ncs{i}),Vn.(ncs{i}),Vc.(ncs{i})]=...
    nuclearSize(1,'static',{ncs{i}, M.(ncs{i})},'mitosis2D');
    else
[An.(ncs{i}),Am.(ncs{i}),Acs.(ncs{i}),Vn.(ncs{i}),Vc.(ncs{i})]=...
    nuclearSize(1,'static',{ncs{i}, M.(ncs{i})},'interphase2D');
    end
end
%% MAIN
if ~isempty(Params)
for rows = 1:size(Params,1)
%     try
%     if rows > 1
%         disp([E(rows-1) penalty(rows-1,:)])
%     end
    Params_ = 10.^(Params(rows,:));
    Params_ = num2cell(Params_);
    [lambdaU,lambdaW,lambdaV,sigmaU,sigmaW,sigmaV,muU,muW,muV,...
        gamma,psi,alpha,beta,phi,beta0,nu,omega,eta,rho,epsilon,beta2,phi2]=Params_{:};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2x Dorsal %%%%%%%%%%%%%%%%%%%%%%%%%%%%


        [x1,x2]=meshgrid(linspace(0,1,N.(nc)),linspace(0,1,M.(nc)));
        toll = getToll(x1,x2,phi,[beta2/beta,phi2]);


    y0(:,:,1) = zeros(M.NC10,N.NC10,1);
    y0(:,:,2) = zeros(M.NC10,N.NC10,1);
    y0(:,:,3) = args.dl*ones(M.NC10,N.NC10,1);
    y0(:,:,4) = args.dl*ones(M.NC10,N.NC10,1);
    y0(:,:,5) = args.cact*(1/alpha)*ones(M.NC10,N.NC10,1);
    y0(:,:,6) = args.cact*(1/alpha)*ones(M.NC10,N.NC10,1);
    y0(:,:,7) = (beta/rho).*Vc.NC10*rho/(Acs.NC10*epsilon*eta)*toll;
    y0(:,:,8) = zeros(M.NC10,N.NC10,1);
    
    
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
%     protein2x=protein;
% c = [protein.dlNuc.NC14(:,:,end); flipud(protein.dlNuc.NC14(:,:,end))];
% c = c(:,round(0.5:0.5:size(c,2)));
% c=c';
% sna = double(c>=1.05);

    if args.plot        
%         [xe,ye,ze] = ellipsoid(0,0,0,1,1,2,length(c)-1);
%         subplot(1,2,1)
%         surf(ze,xe,ye,c)
%         view([0 0])
%         title('free dl')
%         shading interp
%         axis equal
% 
%         subplot(1,2,2)
%         surf(ze,xe,ye,sna)
%         view([0 0])
%         title('sna exp')
%         shading interp
%         axis equal
%         
        
        figure; subplot(1,2,1)
        cmax = max(protein.nuclearDorsal.NC14(:));
%         for j = 1:2:length(ncs)
        for j = 1:2:length(ncs)
            d = protein.dlNuc.(ncs{j});
            [xe,ye]=meshgrid(linspace(0,1,N.(ncs{j})),linspace(0,1,M.(ncs{j})));
            for tim = 1:size(d,3)
                c = d(:,:,tim);
                surf(xe,ye,c)
                title('free dl')
                shading flat
                caxis([0 cmax])
                zlim([0 cmax])
                drawnow
            end
        end
        subplot(1,2,2)
        sna = double(c>=sna_threshold);
        surf(xe,ye,sna)
        shading flat
        view([0 90])
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1x Dorsal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     [x1,x2]=meshgrid(linspace(0,1,M.NC10),linspace(0,1,N.NC10));
%         toll = getToll(x1,x2,phi);
%     y0 = zeros(M.NC10,N.NC10,8);
%     y0(:,:,1) = zeros(M.NC10,N.NC10,1);
%     y0(:,:,2) = zeros(M.NC10,N.NC10,1);
%     y0(:,:,3) = args.dl*ones(M.NC10,N.NC10,1);
%     y0(:,:,4) = args.dl*ones(M.NC10,N.NC10,1);
%     y0(:,:,5) = args.cact*(1/alpha)*ones(M.NC10,N.NC10,1);
%     y0(:,:,6) = args.cact*(1/alpha)*ones(M.NC10,N.NC10,1);
%     y0(:,:,7) = (beta/rho).*Vc.NC10*rho/(Acs.NC10*epsilon*eta)*toll;
%     y0(:,:,8) = zeros(M.NC10,N.NC10,1);
%     
%     nt=1;
%     for i = 1:9                                                         
%         if mod(i,2)==0                                                  
%             nc = ncs{i};                                                
%             mitosis;                                                    
%         else                                                            
%             nc = ncs{i};                                                
%             interphase;                                                 
%         end                                                             
%     end                                                                 
%         
%     for j=1:length(ncs)
%         protein.('nuclearDorsal').(ncs{j})=protein.dlNuc.(ncs{j})+...
%             protein.dlCactNuc.(ncs{j});
%         protein.('totalDorsal').(ncs{j})...
%             = (Vn.(ncs{j})*(protein.dlNuc.(ncs{j})+...
%             protein.dlCactNuc.(ncs{j}))+Vc.(ncs{j})*(protein.dlCyt.(ncs{j})+...
%             protein.dlCactCyt.(ncs{j})))/(Vc.(ncs{j})+Vn.(ncs{j}));
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
D = [reshape(protein.nuclearDorsal.NC11(:,end,:),[],1);
    reshape(protein.nuclearDorsal.NC12(:,end,:),[],1);
    reshape(protein.nuclearDorsal.NC13(:,end,:),[],1);
    reshape(protein.nuclearDorsal.NC14(:,end,:),[],1)];
    Beta = mean(D.*C)/mean(D.^2);
    E(rows) = sqrt(sum(((C-Beta*D)./dC).^2)/length(D));
%     snaE = zeros(50);
%     snaE(1:12,:) = 1;
%     snaE(18:50,[1 2 3 48 49 50])=1;
%     E = sum(sum((snaE-sna).^2));
end
end
%% Output

pnames = {'lambdaU','lambdaW','lambdaV','sigmaU','sigmaW','sigmaV',...
    'muU','muW','muV','gamma','psi','alpha','beta','phi','beta0','nu',...
    'omega','eta','rho','epsilon','beta2','phi2'};
argsout = {E,penalty,protein,time,Beta,pnames};

for i = 1:nargout
    varargout{i} = argsout{i};
end

%% nested functions
    function interphase
        %update variables
        [x1,x2]=meshgrid(linspace(0,1,N.(nc)),linspace(0,1,M.(nc)));
        toll = getToll(x1,x2,phi,[beta2/beta,phi2]);
        
        h=1;
        g=1;
        e1 = ones(1,M.(nc));
        e2 = ones(1,N.(nc));
        
        e1(end) = 0; e1(end-1) = 2; e1 = repmat(e1,N.(nc),1); e1 = e1(:);
        e2(end) = 0; e2(end-1) = 2; e2 = repmat(e2',M.(nc),1); e2 = e2(:);
        e3 = ones(M.(nc)*N.(nc),1);
        e4 = flipud(e2);
        e5 = flipud(e1);
        P = spdiags([e1/g^2 e2/h^2 -2*e3/g^2-2*e3/h^2 e4/h^2 e5/g^2],...
            [-N.(nc) -1 0 1 N.(nc)],M.(nc)*N.(nc),M.(nc)*N.(nc));
        

        options = odeset('RelTol',1e-1,'AbsTol',1e-2,'JPattern',jacpat(M.(nc),N.(nc)));
        

        
        % interpolating (nuc = cyt)
        [a,b,~]=size(y0);
        [x1o,x2o] = meshgrid(linspace(0,1,b),linspace(0,1,a));
        
        y0_ = y0;
        y0 = zeros(M.(nc),N.(nc),8);
        
        for n = 1:2:6
            y0(:,:,n) = interp2(x1o,x2o,y0_(:,:,n+1),x1,x2);
            y0(:,:,n+1) = interp2(x1o,x2o,y0_(:,:,n+1),x1,x2);
        end
        for n = 7:8
            y0(:,:,n)=interp2(x1o,x2o,y0_(:,:,n),x1,x2);
        end
        y0(isnan(y0))=0;
        y0_ = zeros(size(y0(:)));
        nn = length(y0(:))/8;
        for n = 1:8
            temp = y0(:,:,n);
            y0_((n-1)*nn+1:(n*nn),1) = temp(:);
        end
        y0 = y0_;
        if times
            [t,y] = ode15s(@int_fun,tspan.(nc),y0,options);
            [time.(xs{nt}),time.(ys{nt}),time.(ts{nt})] = meshgrid(linspace(0,1,M.(nc)),linspace(0,1,N.(nc)),t);
            nt=nt+1;
        else
            [~,y] = ode15s(@int_fun,tspan.(nc),y0,options);
        end
        
        y = y';
        
        for i_ = 1:length(names)
            temp = y((i_-1)*M.(nc)*N.(nc)+1:(i_)*M.(nc)*N.(nc),:);
            tempNew = zeros(M.(nc),N.(nc),size(temp,2));
            for j_ = 1:size(temp,2)
                tempNew(:,:,j_) = reshape(temp(:,j_),M.(nc),N.(nc));
            end
            protein.(names{i_}).(nc)=tempNew;
        end
        
        % update y0 for mitosis (mixing)
        y0_ = y(:,end);
        y0 = zeros(M.(nc),N.(nc),8);
        for i_ =[6 5 4 3 2 1]
            if mod(i_,2)==0
                y0(:,:,i_) = reshape((Vn.(nc)*y0_((i_-2)*M.(nc)*N.(nc)+1:...
                    (i_-1)*M.(nc)*N.(nc)) + ...
                    Vc.(nc)*y0_((i_-1)*M.(nc)*N.(nc)+1:...
                    (i_)*M.(nc)*N.(nc)))/(Vn.(nc)+Vc.(nc)),[M.(nc),N.(nc),1]);
            end
            
        end
        for i_ = [7 8]
                y0(:,:,i_) = reshape((Vn.(nc)*y0_((i_-2)*M.(nc)*N.(nc)+1:...
                    (i_-1)*M.(nc)*N.(nc)) + ...
                    Vc.(nc)*y0_((i_-1)*M.(nc)*N.(nc)+1:...
                    (i_)*M.(nc)*N.(nc)))/(Vn.(nc)+Vc.(nc)),[M.(nc),N.(nc),1]);
        end
        
        y0_ = zeros(size(y0(:)));
        nn = length(y0(:))/8;
        for n = 1:8
            temp = y0(:,:,n);
            y0_((n-1)*nn+1:(n*nn),1) = temp(:);
        end
        y0 = y0_;
        
    end
    function F = int_fun(~,y)
        
        un = y(1:M.(nc)*N.(nc));
        uc = y(M.(nc)*N.(nc)+1:2*M.(nc)*N.(nc));
        wn = y(2*M.(nc)*N.(nc)+1:3*M.(nc)*N.(nc));
        wc = y(3*M.(nc)*N.(nc)+1:4*M.(nc)*N.(nc));
        vn = y(4*M.(nc)*N.(nc)+1:5*M.(nc)*N.(nc));
        vc = y(5*M.(nc)*N.(nc)+1:6*M.(nc)*N.(nc));
        yy = y(6*M.(nc)*N.(nc)+1:7*M.(nc)*N.(nc));
        xx = y(7*M.(nc)*N.(nc)+1:8*M.(nc)*N.(nc));
        
        f1 = (An.(nc)*(sigmaU*uc-muU*un)+Vn.(nc)*(beta0*wn...
            -gamma*un.*vn))/(Vn.(nc));
        f2 = (lambdaU*Am.(nc)*P*uc+omega*epsilon*Acs.(nc)*xx...
            +Vc.(nc)*(beta0*wc-gamma*uc.*vc)...
            -An.(nc)*(sigmaU*uc-muU*un))/(Vc.(nc));
        f3 = (An.(nc)*(sigmaW*wc-muW*wn)-Vn.(nc)*(beta0*wn...
            -gamma*un.*vn))/(Vn.(nc));
        f4 = (lambdaW*Am.(nc)*P*wc-Acs.(nc)*epsilon*(eta*wc.*yy-nu*xx)...
            -Vc.(nc)*(beta0*wc-gamma*uc.*vc)...
            -An.(nc)*(sigmaW*wc-muW*wn))/(Vc.(nc));
        f5 = (An.(nc)*(sigmaV*vc-muV*vn)+Vn.(nc)*psi*(beta0*wn...
            -gamma*un.*vn))/(Vn.(nc));
        f6 = (lambdaV*Am.(nc)*P*vc+psi*(omega*epsilon*Acs.(nc)*xx...
            +Vc.(nc)*(beta0*wc-gamma*uc.*vc))+1-alpha*Vc.(nc)*vc...
            -An.(nc)*(sigmaV*vc-muV*vn))/(Vc.(nc));
        f7 = (nu+omega)*xx-eta*yy.*wc+beta*toll(:)-rho*yy;
        f8 = eta*yy.*wc-(nu+omega)*xx;
        
        F = [f1;f2;f3;f4;f5;f6;f7;f8];
    end

    function mitosis
        % update variables
        if times
            [t,y] = ode15s(@mit_fun,tspan.(nc),y0,options);
            [time.(xs{nt}),time.(ys{nt}),time.(ts{nt})]...
                = meshgrid(linspace(0,1,M.(nc)),linspace(0,1,N.(nc)),t);
            nt=nt+1;
        else
            [~,y] = ode15s(@mit_fun,tspan.(nc),y0,options);
        end
        y = y';
        y0 = zeros(M.(nc),N.(nc),8);
        for i_ = 1:length(names)
            temp = y((i_-1)*M.(nc)*N.(nc)+1:(i_)*M.(nc)*N.(nc),:);
            tempNew = zeros(M.(nc),N.(nc),size(temp,2));
            for j_ = 1:size(temp,2)
                tempNew(:,:,j_) = reshape(temp(:,j_),M.(nc),N.(nc));
            end
            protein.(names{i_}).(nc)=tempNew;
            
            % update y0 for interphase (interpolation)
            y0(:,:,i_) = tempNew(:,:,end);
        end
    end

    function F = mit_fun(~,y)
        
        uc = y(M.(nc)*N.(nc)+1:2*M.(nc)*N.(nc));
        wc = y(3*M.(nc)*N.(nc)+1:4*M.(nc)*N.(nc));
        vc = y(5*M.(nc)*N.(nc)+1:6*M.(nc)*N.(nc));
        yy = y(6*M.(nc)*N.(nc)+1:7*M.(nc)*N.(nc));
        xx = y(7*M.(nc)*N.(nc)+1:8*M.(nc)*N.(nc));
        
        
        f1 = zeros(M.(nc)*N.(nc),1);
        f2 = (lambdaU*Am.(nc)*P*uc+omega*epsilon*Acs.(nc)*xx...
            -gamma*Vc.(nc)*uc.*vc+beta0*Vc.(nc)*wc)/(Vc.(nc));
        f3 = zeros(M.(nc)*N.(nc),1);
        f4 = (lambdaW*Am.(nc)*P*wc-eta*epsilon*Acs.(nc)*wc.*yy...
            +nu*epsilon*Acs.(nc)*xx+gamma*Vc.(nc)*uc.*vc...
            -beta0*Vc.(nc)*wc)/(Vc.(nc));
        f5 = zeros(M.(nc)*N.(nc),1);
        f6 = (lambdaV*Am.(nc)*P*vc+psi*(omega*epsilon*Acs.(nc)*xx...
            -gamma*Vc.(nc)*uc.*vc+beta0*Vc.(nc)*wc)+1-...
            alpha*Vc.(nc)*vc)/(Vc.(nc));
        f7 = (nu+omega)*xx-eta*yy.*wc+beta*toll(:)-rho*yy;
        f8 = eta*yy.*wc-(nu+omega)*xx;
        
        F = [f1;f2;f3;f4;f5;f6;f7;f8];
        
    end
end

function pat = jacpat(M,N)
e1 = ones(1,M); e1 = repmat(e1,N,1); e1 = e1(:);
e2 = ones(1,N); e2 = repmat(e2',M,1); e2 = e2(:);
e3 = ones(M*N,1);
e4 = flipud(e2);
e5 = flipud(e1);
P = spdiags([e1 e2 e3 e4 e5],[-N -1 0 1 N],M*N,M*N);
D = speye(M*N);
Z = 0.*D;

pat = [D D D Z D Z Z Z
    D P Z D Z D Z D 
    D Z D D D Z Z Z
    Z D D P Z D D D
    D Z D Z D D Z Z
    Z D Z D D P Z D
    Z Z Z D Z Z D D
    Z Z Z D Z Z D D];
end

function toll = getToll(x,y,phi,vars)
if length(vars) == 2
    toll = vars(1)*exp(-0.5*((y)./(vars(2))).^2)+exp(-0.5*((x)./(phi)).^2);
else
    error('getToll variable "vars" needs to have 2 entries')
end
end