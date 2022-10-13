% script_ambrosi2014_shuttling_Diff
%
% This script looks at the effect of diffusion just like we did for the
% analysis of Fig 5C of Ambrosi, but using the shuttling model with our
% best parameter sets.
%

close all

addpath('C:\Users\gtreeves\Documents\Dropbox\Matlab\Dorsal\MDO\shuttling')
addpath('C:\Users\gtreeves\Documents\Dropbox\Matlab\Dorsal\MDO\shuttling\Other Files')
load('results999999_dsat_jp2_mdo.mat','xb')

An = 1; M=51; Vn = 1; Vc = 2.7665;
ps2ch = { [1 2 3], [1 2], [ 1 3], 1,[2 3],  2, 3};
titles = {'Global','dl + dlCact','dl + cact',  'dl','dlCact + cact', 'dlCact', 'Cact'};
x = linspace_offset(0,1,51);
x=x';
for k = 1:5
P = xb(k,:);
figure1 = figure; set(gcf,'Name','Effects of diffusion')
% h2 = figure; set(gcf,'Name','Corresponding Fluxes')

nn = 8;
rgb = colormap('cool');
rgb = rgb(round(linspace(20,length(rgb),nn))',:);
for i = 1:length(ps2ch)

    Pd = P;
    
    figure(figure1)
    for j = 1:nn
        Pd(ps2ch{i}) = P(ps2ch{i})-(j-1)*0.25;
        
        subplot(4,7,i)
        hold on
        [~,~,PROTEIN,~,~,~] = dsat_jp2(Pd);
        
        dn1 = (PROTEIN.dlNuc.NC14(:,end));
        dn1 = dn1-min(dn1(:)); dn1 = dn1/max(dn1(:));
        plot(x,dn1,'Color',rgb(j,:))
        title(titles{i})
        
        %         figure(h2)
        subplot(4,7,i+7)
        hold on
        dC = -diff((PROTEIN.dlCyt.NC14(:,end)))*Pd(1);
        plot(x(1:end-1),dC,'Color',rgb(j,:))
        title('dl flux')
        
        subplot(4,7,i+14)
        hold on
        dcC = -diff((PROTEIN.dlCactCyt.NC14(:,end)))*Pd(2);
        plot(x(1:end-1),dcC,'Color',rgb(j,:))
        title('dlCact flux')
        subplot(4,7,i+21)
        hold on
        cC = -diff((PROTEIN.cactCyt.NC14(:,end)))*Pd(3);
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

end