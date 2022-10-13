% Ambrosi2014 - Kanodia2009
%
% This script looks at the effect of changing both Toll strength and
% k_import, or both of those with k_export too.
%
%
close all
% From Table 1, Ambrosi2014
L = 241.5; %um
Er = 102.4; %um
Eh = 25; %um
Tn = 6000; 
R=15000;
% R = 150000;
S = 4000;
xi = 3;

Lambda=4;
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


%% Kanodia model

An14 = 160.1618;
Vn14 = 186.034;
Am14 = 264.7;
T = 1;

% f1 = figure;
% hold on
% set(gcf,'Position',[1 5 1280 633])
% set(gcf,'Name','Kanodia Model')

x = linspace_offset(0,1,51);
x=x';
x = [flipud(-x);x];
figure; set(gcf,'Name','Large perturbations of ki')
for i = 1:5
% for i = 3
del = 10^(i);
    ki = 4*del;
     
lambdaU = Lambda*Am14*T/Vn14;
lambdaW = Lambda*Am14*T/Vn14;
lambdaV = Lambda*Am14*T/Vn14;
sigmaU = ki*An14*T/Vn14;
muU = ke*An14*T/Vn14;
gamma = kb*Cact0*T;
psi = Cdl0/Cact0;
alpha = kDeg*T;

beta = R*(i/5)*T/L^xi;
phi = S/L^xi;

         
P = log10([lambdaU,lambdaW,lambdaV,sigmaU,muU,...
            gamma,psi,alpha,beta,phi,xi]);
An = 1; 
M=51; Vn = 1; Vc = 2.7665;
[ERROR,PENALTY,PROTEIN,TIME,BETA] = kanodia_log_nest_JPattern(P);
dlNuc = PROTEIN.dlNuc.NC14;
dlCyt = PROTEIN.dlCyt.NC14;
dlCact = PROTEIN.dlCactCyt.NC14;
Cact = PROTEIN.cactCyt.NC14;

 nc = {'NC10','NC11','NC12','NC13','NC14'};
for k = 1:5
subplot(5,5,(i-1)*5+k), hold on
dn = (PROTEIN.dlNuc.(nc{k}));
for j = 1:7
plot(dn(j,:)','Color',[0 j-1 j-1]/6)
end
title(nc{k})
end

end

figure; set(gcf,'Name','Realistic perturbations of ki, ke and lambda')
for i = 1:5
% for i = 3
del = 10^((3-i));
    ki = 4*del;
    ke = 1*del;
    Lambda = 3*del;
    
lambdaU = Lambda*Am14*T/Vn14;
lambdaW = Lambda*Am14*T/Vn14;
lambdaV = Lambda*Am14*T/Vn14;
sigmaU = ki*An14*T/Vn14;
muU = ke*An14*T/Vn14;
gamma = kb*Cact0*T;
psi = Cdl0/Cact0;
alpha = kDeg*T;

beta = R*(i/5)*T/L^xi;
phi = S/L^xi;

         
P = log10([lambdaU,lambdaW,lambdaV,sigmaU,muU,...
            gamma,psi,alpha,beta,phi,xi]);
An = 1; 
M=51; Vn = 1; Vc = 2.7665;
[ERROR,PENALTY,PROTEIN,TIME,BETA] = kanodia_log_nest_JPattern(P);
dlNuc = PROTEIN.dlNuc.NC14;
dlCyt = PROTEIN.dlCyt.NC14;
dlCact = PROTEIN.dlCactCyt.NC14;
Cact = PROTEIN.cactCyt.NC14;

 nc = {'NC10','NC11','NC12','NC13','NC14'};
for k = 1:5
subplot(5,5,(i-1)*5+k), hold on
dn = (PROTEIN.dlNuc.(nc{k}));
for j = 1:7
plot(dn(j,:)','Color',[0 j-1 j-1]/6)
end
title(nc{k})
end

end
% formatting
% set(gca,'View',[0 0])

% figure(f2)
% set(gca,'YScale','log')
% set(gca,'ZScale','lin')

%% ke
% figure
% set(gcf,'Position',[1 5 1280 633])
% set(gcf,'Name','Kanodia Model')
% subplot(1,2,1)
% hold on
% title('nuclear dl')
% subplot(1,2,2)
% hold on
% title('dl/Cact')
% 
% for i = 1:6
%     ke = 10^(i-2);
% lambdaU = LambdaU*Am14*T/Vn14;
% lambdaW = LambdaW*Am14*T/Vn14;
% lambdaV = LambdaV*Am14*T/Vn14;
% sigmaU = ki*An14*T/Vn14;
% muU = ke*An14*T/Vn14;
% gamma = kb*Cact0*T;
% psi = Cdl0/Cact0;
% alpha = kDeg*T;
% 
% beta = R*T/L^xi;
% phi = S/L^xi;
% 
%          
% P = log10([lambdaU,lambdaW,lambdaV,sigmaU,muU,...
%             gamma,psi,alpha,beta,phi,xi]);
% 
% 
% [ERROR,PENALTY,PROTEIN,TIME,BETA] = kanodia_log_nest_JPattern(P);
% 
% subplot(1,2,1)
% plot3(linspace_offset(-1,1,nn),10^(i-2)*ones(1,nn),[flipud(PROTEIN.dlNuc.NC14(:,end));PROTEIN.dlNuc.NC14(:,end)])
% set(gca,'YScale','log')
% subplot(1,2,2)
% plot3(linspace_offset(-1,1,nn),10^(i-2)*ones(1,nn),[flipud(dlCact(:,end));dlCact(:,end)])
% set(gca,'YScale','log')
% 
% end
% subplot(1,2,1)
% xlabel('DV Coordinate','FontSize',18)
% ylabel('fold change, k_{e}','FontSize',18)
% zlabel('NC 14 (AU)')
% view(-30,30)
% 
% subplot(1,2,2)
% xlabel('DV Coordinate','FontSize',18)
% ylabel('fold change, k_{e}','FontSize',18)
% zlabel('NC 14 (AU)')
% view(-30,30)