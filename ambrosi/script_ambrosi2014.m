% Ambrosi2014 - Kanodia2009
%
%
% This script looks at the basic output of the ambrosi model.


addpath('C:\Users\gtreeves\Documents\Dropbox\Matlab\Dorsal\MDO\deconvolution')
addpath('C:\Users\gtreeves\Documents\Dropbox\Matlab\Dorsal\MDO\deconvolution\Other Files')

% From Table 1, Ambrosi2014
L = 241.5; %um
Er = 102.4; %um
Eh = 25; %um
Tn = 6000; 
R=15000;
S = 4000;
xi = 2.5;
LambdaU=4; LambdaW=4; LambdaV=4;
ki =4;
ke =1;
Pcact =50;
kDeg=1;
kb=0.02;
Cdl0 = 36;
Cdlcact0 = 30;
Cact0 = 36;

r = 3.08; %um
n = 92/2;
t = 65; % min

An = 4*pi*r*r;
Am = Eh*Eh*2;
Vn = 4/3*pi*r*r*r;
Vc = L*Eh*Eh/n-Vn;


% [ y ] = fun_abrosi(t,n,L,An,ki,ke,Am,Vn,LambdaU,Vc,kb,kDeg,Pcact,R,S,xi);
% y=y';
% dln = y(1:n,:);
% dlc = y(n+1:2*n,:);
% dlcact = y(2*n+1:3*n,:);
% cact = y(3*n+1:4*n,:); 
% 
% figure
% set(gcf,'Position',[1 711 1280 634])
% set(gcf,'Name','Ambrosi Model')
% subplot(2,3,1)
% plotyy(linspace(-1,1,n*2),[flipud(dln(:,end));dln(:,end)],...
%     linspace(-1,1,n*2),[flipud(dlcact(:,end));dlcact(:,end)])
% legend('dl n','dlCact cyt')
% 
% subplot(2,3,2)
% surf([flipud(dln); dln] )
% shading interp
% title('dl nuc')
% 
% subplot(2,3,3)
% surf([flipud(dlc); dlc] )
% shading interp
% title('dl cyt')
% 
% subplot(2,3,5)
% surf([flipud(dlcact); dlcact] )
% shading interp
% title('dl/Cact')
% 
% subplot(2,3,6)
% surf([flipud(cact); cact] )
% shading interp
% title('cact')

%% Kanodia model

An14 = 160.1618;
Vn14 = 186.034;
Am14 = 264.7;
T = 1;


lambdaU = LambdaU*Am14*T/Vn14;
lambdaW = LambdaW*Am14*T/Vn14;
lambdaV = LambdaV*Am14*T/Vn14;
sigmaU = ki*An14*T/Vn14;
muU = ke*An14*T/Vn14;
gamma = kb*Cact0*T;
psi = Cdl0/Cact0;
alpha = kDeg*T;

beta = R*T/L^xi;
phi = S/L^xi;

         
P = log10([lambdaU,lambdaW,lambdaV,sigmaU,muU,...
            gamma,psi,alpha,beta,phi,xi]);


[ERROR,PENALTY,PROTEIN,TIME,BETA] = kanodia_log_nest_JPattern(P);

figure
set(gcf,'Position',[1 5 1280 633])
set(gcf,'Name','Kanodia Model')
subplot(2,3,1)
plotyy(linspace(-1,1,102),[flipud(PROTEIN.dlNuc.NC14(:,end));PROTEIN.dlNuc.NC14(:,end)],...
    linspace(-1,1,102),[flipud(PROTEIN.dlCactCyt.NC14(:,end));PROTEIN.dlCactCyt.NC14(:,end)])
legend('dl n','dlCact cyt')


subplot(2,3,2)
surf([flipud(PROTEIN.dlNuc.NC14); PROTEIN.dlNuc.NC14] )
shading interp
title('dl nuc')

subplot(2,3,3)
surf([flipud(PROTEIN.dlCyt.NC14); PROTEIN.dlCyt.NC14] )
shading interp
title('dl cyt')

subplot(2,3,5)
surf([flipud(PROTEIN.dlCactCyt.NC14); PROTEIN.dlCactCyt.NC14] )
shading interp
title('dl/Cact')

subplot(2,3,6)
surf([flipud(PROTEIN.cactCyt.NC14); PROTEIN.cactCyt.NC14] )
shading interp
title('cact')
