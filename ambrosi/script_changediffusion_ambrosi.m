% Ambrosi2014 - Kanodia2009
%
% Formerly titled script_ambrosi2014_Diff
%
% This script evaluates the effect of changing certain sets of diffusion
% parameters on the Ambrosi model.  Specifically, we are looking at Figure
% 5C from Ambrosi, which examined the gyn mutant.  In that subfigure, it
% appears that changing the value of Gamma, which controls the
% intercompartmental flux of all species, from 0.03 to 2 causes the
% gradient to become sharper and taller.  This would argue for some
% shuttling-like effect, or that a model that did not take shuttling into
% consideration could still explain the diffusion results we are giving.
%
% The problem is that, even if a model does not "consider shuttling", it
% still may have the shuttling phenomenon.  So in this script, we show that
% the shuttling phenomenon is there.  If we systematically decrease certain
% groups of diffusion parameters, we find that as long as dl diffusion is
% decreased, the gradient sharpens (as you'd expect), even if you also
% decrease dl/Cact complex diffusion in the same way.  But if you decrease
% dl/Cact complex diffusion without changing dl diffusion, then the
% gradient widens (as you'd expect from a shuttling model).  To show that
% shuttling is happening, if we examine the fluxes under these conditions,
% the dl/Cact flux drops precipitously, without much drop in the dl flux.
% This means the gradient expansion comes from the lack of shuttling.

% load Mat/gyn
load Mat/Ambrosi2014_baseparams
load Mat/Ambrosi2014_3rdcol

ps2ch = { [1 2 3], [1 2], [ 1 3], 1,[2 3],  2, 3};
titles = {'Global','dl + dlCact','dl + cact',  'dl','dlCact + cact', 'dlCact', 'Cact'};
% x = linspace_offset(0,1,51);
% x=x';

figure1 = figure; set(gcf,'Name','Effects of diffusion')
% h2 = figure; set(gcf,'Name','Corresponding Fluxes')
Lambda_orig = 2*ones(3,1);

nn = 9;
rgb = colormap('cool');
rgb = rgb(round(linspace(20,length(rgb),nn))',:);
for i = 1:length(ps2ch)

    
    
    figure(figure1)
    for j = 1:nn
        Lambda = Lambda_orig;
        Lambda(ps2ch{i}) = 10.^(log10(Lambda(ps2ch{i}))-(j-1)*0.25);
        
        subplot(4,7,i)
        hold on
% disp(Gamma)
        [ dn1,dc1,dcc1,cc1 ] = ambrosifun( El,Er,Eh,Tn,R,S,xi,Lambda,ki,ke,Pcact,kDeg,kb,r,n,tspan,Cdl0,Cdlcact0,Cact0);

        dn1 = dn1(n/2+1:end,end);
        dc1 = dc1(n/2+1:end,end);
        dcc1 = dcc1(n/2+1:end,end);
        cc1 = cc1(n/2+1:end,end);
        
        dn1 = dn1-min(dn1(:)); dn1 = dn1/max(dn1(:));
        plot(1:length(dn1),dn1,'Color',rgb(j,:))
        title(titles{i})
        
        
        subplot(4,7,i+7)
        hold on
        dC = -diff((dc1(:,end)))*Lambda(1);
        plot(1:1:(length(dc1)-1),dC,'Color',rgb(j,:))
        title('dl flux')
        
        subplot(4,7,i+14)
        hold on
        dcC = -diff((dcc1(:,end)))*Lambda(2);
        plot(1:1:(length(dc1)-1),dcC,'Color',rgb(j,:))
        title('dlCact flux')
        subplot(4,7,i+21)
        hold on
        cC = -diff((cc1(:,end)))*Lambda(3);
        plot(1:1:(length(dc1)-1),cC,'Color',rgb(j,:))
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

