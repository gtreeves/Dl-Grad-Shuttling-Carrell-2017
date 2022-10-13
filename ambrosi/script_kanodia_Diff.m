% Ambrosi2014 - Kanodia2009
% close all
% From Table 1, Ambrosi2014


addpath('C:\Users\gtreeves\Documents\Dropbox\Matlab\Dorsal\MDO\deconvolution')
addpath('C:\Users\gtreeves\Documents\Dropbox\Matlab\Dorsal\MDO\deconvolution\Other Files')


L = 241.5; %um
Er = 102.4; %um
Eh = 25; %um
Tn = 6000;
R=15000;
% R = 150000;
S = 4000;
xi = 3;

Lambda=2;
ki =40/100;
ke =.44/100;
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
An = 1; M=51; Vn = 1; Vc = 2.7665;
ps2ch = { [1 2 3], [1 2], [ 1 3], 1,[2 3],  2, 3};
titles = {'Global','dl + dlCact','dl + cact',  'dl','dlCact + cact', 'dlCact', 'Cact'};
x = linspace_offset(0,1,51);
x=x';

figure1 = figure; set(gcf,'Name','Effects of diffusion')
% h2 = figure; set(gcf,'Name','Corresponding Fluxes')

nn = 9;
rgb = colormap('cool');
rgb = rgb(round(linspace(20,length(rgb),nn))',:);
for i = 1:length(ps2ch)
    % for i = 3
    
    lambdaU = Lambda*Am14*T/Vn14;
    lambdaW = Lambda*Am14*T/Vn14;
    lambdaV = Lambda*Am14*T/Vn14;
    sigmaU = ki*An14*T/Vn14;
    muU = ke*An14*T/Vn14;
    gamma = kb*Cact0*T;
    psi = Cdl0/Cact0;
    alpha = kDeg*T;
    
    beta = R*T/L^xi;
    phi = S/L^xi;
    
    
    P = log10([lambdaU,lambdaW,lambdaV,sigmaU,muU,...
        gamma,psi,alpha,beta,phi,xi]);
    Pd = P;
    
    figure(figure1)
    for j = 1:nn
        Pd(ps2ch{i}) = P(ps2ch{i})-(j-1)*0.25;
        
        subplot(4,7,i)
        hold on
        [~,~,PROTEIN,~,~,~] = kanodia_log_nest_JPattern(Pd);
        
%         dn1{j} = (PROTEIN.dlNuc.NC14(:,end));
        dn1 = (PROTEIN.dlNuc.NC14(:,end));
        dn1 = dn1-min(dn1(:)); dn1 = dn1/max(dn1(:));
        plot(x,dn1,'Color',rgb(j,:))
        title(titles{i})
        
        %         figure(h2)
        subplot(4,7,i+7)
        hold on
        dC = -diff((PROTEIN.dlCyt.NC14(:,end)))*10.^Pd(1);
        plot(x(1:end-1),dC,'Color',rgb(j,:))
        title('dl flux')
        
        subplot(4,7,i+14)
        hold on
        dcC = -diff((PROTEIN.dlCactCyt.NC14(:,end)))*10.^Pd(2);
        plot(x(1:end-1),dcC,'Color',rgb(j,:))
        title('dlCact flux')
        subplot(4,7,i+21)
        hold on
        cC = -diff((PROTEIN.cactCyt.NC14(:,end)))*10.^Pd(3);
        plot(x(1:end-1),cC,'Color',rgb(j,:))
        title('Cact flux')
    end
    
    
end
% Create arrow
annotation(figure1,'arrow',[0.198046875 0.15859375],...
    [0.840791044776119 0.841044776119403]);

% Create arrow
annotation(figure1,'arrow',[0.317187500000001 0.277734375000001],...
    [0.841537313432836 0.841791044776119]);

% Create arrow
annotation(figure1,'arrow',[0.428125000000001 0.388671875000001],...
    [0.840044776119403 0.840298507462687]);

% Create arrow
annotation(figure1,'arrow',[0.542968750000001 0.503515625000001],...
    [0.840044776119403 0.840298507462687]);

% Create arrow
annotation(figure1,'arrow',[0.619140625 0.655859375],...
    [0.83855223880597 0.838059701492537]);

% Create arrow
annotation(figure1,'arrow',[0.733984375000025 0.770703125000025],...
    [0.83855223880597 0.838059701492537]);

% Create textbox
annotation(figure1,'textbox',...
    [0.0500000000000001 0.38955223880597 0.0578125 0.0333283582089555],...
    'String',{'-\lambdadC/dx'},...
    'LineStyle','none',...
    'FontSize',36,...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.0523437500000001 0.615671641791045 0.0578125 0.0333283582089554],...
    'String',{'-\lambdadC/dx'},...
    'LineStyle','none',...
    'FontSize',36,...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.0507812500000001 0.168656716417911 0.0578125 0.0333283582089554],...
    'String',{'-\lambdadC/dx'},...
    'LineStyle','none',...
    'FontSize',36,...
    'FitBoxToText','off');

saveas(gcf,'Mat/script_kanodia_Diff.fig')