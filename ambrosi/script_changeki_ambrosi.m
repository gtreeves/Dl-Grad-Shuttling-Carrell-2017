% Ambrosi2014 - Kanodia2009
%
% This script runs a for loop to change k_import.
%
% close all

addpath('C:\Users\gtreeves\Documents\Dropbox\Matlab\Dorsal\MDO\deconvolution')
addpath('C:\Users\gtreeves\Documents\Dropbox\Matlab\Dorsal\MDO\deconvolution\Other Files')

% From Table 1, Ambrosi2014
L = 241.5; %um
Er = 102.4; %um
Eh = 25; %um
Tn = 6000; 
R=15000;
% R = 150000;
S = 4000;
xi = 3;

LambdaU=0; LambdaW=LambdaU; LambdaV=LambdaU;
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

f1 = figure;
set(gcf,'Position',[1 5 1280 633])
set(gcf,'Name','Kanodia Model')
subplot(2,2,1)
hold on
title('nuclear dl')
subplot(2,2,2)
hold on
title('dl/Cact')
subplot(2,2,3)
hold on
title('dl Cyt')
subplot(2,2,4)
hold on
title('Cact')
f2 = figure;
hold on
x = linspace_offset(0,1,51);
x=x';
x = [flipud(-x);x];
for i = 1:6
% for i = 3
    ki = 4*10^(i-3);
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
An = 1; 
M=51; Vn = 1; Vc = 2.7665;
[ERROR,PENALTY,PROTEIN,TIME,BETA] = kanodia_log_nest_JPattern(P);
SOL= ode15s(@int_fun_no_diff,[0 50000],[zeros(M*2,1); ones(M*2,1)],[],M,sigmaU,muU,...
            gamma,psi,alpha,beta,phi,xi,Vn,Vc,An);
% dlNuc = SOL.y(1:M,:);
% dlCyt = SOL.y(M+1:2*M,:);
% dlCact = SOL.y(2*M+1:3*M,:);
% Cact = SOL.y(3*M+1:4*M,:);

dlNuc = PROTEIN.dlNuc.NC14;
dlCyt = PROTEIN.dlCyt.NC14;
dlCact = PROTEIN.dlCactCyt.NC14;
Cact = PROTEIN.cactCyt.NC14;
figure;
surf(dlNuc), shading interp
figure(f1)
nn = length([flipud(dlNuc(:,end));dlNuc(:,end)]);
subplot(2,2,1)
plot3(linspace_offset(-1,1,nn),10^(i-2)*ones(1,nn),...
    [flipud(dlNuc(:,end));dlNuc(:,end)])

subplot(2,2,2)
plot3(linspace_offset(-1,1,nn),10^(i-2)*ones(1,nn),...
    [flipud(dlCact(:,end));dlCact(:,end)])

subplot(2,2,3)
plot3(linspace_offset(-1,1,nn),10^(i-2)*ones(1,nn),...
    [flipud(dlCyt(:,end));dlCyt(:,end)])

subplot(2,2,4)
plot3(linspace_offset(-1,1,nn),10^(i-2)*ones(1,nn),...
    [flipud(Cact(:,end));Cact(:,end)])


figure(f2)
dlNuc = [flipud(dlNuc(:,end));dlNuc(:,end)];
dlCact=[flipud(dlCact(:,end));dlCact(:,end)];
A = kb*Pcact*ke/kDeg/ki*dlNuc; A=A-min(A(:)); A = A/max(A(:));
y = beta./(phi+abs(x).^xi);
B = dlCact.*y; B= B-min(B(:)); B=B/max(B(:));
C = 1./(1+muU/sigmaU+gamma*1*muU/alpha/sigmaU/beta*(phi+abs(x).^xi));
C= C-min(C(:)); C=C/max(C(:));
plot3(x,10^(i-2)*ones(1,nn),A,':b','LineWidth',3)
plot3(x,10^(i-2)*ones(1,nn),B,'r')
plot3(x,10^(i-2)*ones(1,nn),C,'G')
legend('dlNuc scaled','Toll Signal * dlCact','Analytical Solution')
end
% formatting
figure(f1)
for i = 1:4
subplot(2,2,i)
set(gca,'YScale','log')
xlabel('DV Coordinate','FontSize',18)
ylabel('fold change, k_{i}','FontSize',18)
zlabel('NC 14 (AU)')
view(-30,30)
end

figure(f2)
set(gca,'YScale','log')
set(gca,'ZScale','lin')

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