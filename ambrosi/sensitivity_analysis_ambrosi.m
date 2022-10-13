%% sensitivity analysis
load('results999999_dsat_jp2_mdo.mat','xb')
close all
xb = xb(1:150,:);
xb_orig = xb;
b = 1;
sens = zeros(length(xb),1); err=sens;
for set = 1:150
    xb = xb_orig(set,:);
    nn = 8;
    ps2ch = {[4,7]};
    del = -0.5;
    % del = repmat(-2:1:2,5,1);
    % figure
    
    n=1;
    for i = length(ps2ch)
        try
            xb_ = xb;
            [e1,~,protein,~,~,~,~] = dsat_jp2(xb_);
            Z1 =  [flipud((protein.nuclearDorsal.NC14(:,end))); (protein.nuclearDorsal.NC14(:,end))];
            [~,~,~,sigma1] = fitgauss_noM(linspace(-1,1,102)',Z1);
            xb_(1,ps2ch{i}) = xb(1,ps2ch{i}) + del;
            [~,~,protein,~,~,~,~] = dsat_jp2(xb_);
            Z2 =  [flipud((protein.nuclearDorsal.NC14(:,end))); (protein.nuclearDorsal.NC14(:,end))];
            [e2,~,~,sigma2] = fitgauss_noM(linspace(-1,1,102)',Z2);
        catch ME
            disp(ME.message)
            sigma1 = NaN; sigma2 = NaN;
        end
%         Z1 = Z1-min(Z1(:)); Z1 = Z1/max(Z1(:));
%         Z2 = Z2-min(Z2(:)); Z2 = Z2/max(Z2(:));
%         sens(set) = sum((Z2-Z1).^2);
        sens(set) = abs((sigma2-sigma1)/sigma1);
        err(set) = abs((e2-e1)/e1);
    end
    
    n=n+1;
    disp(set)
end

%% Analyze the 10% most sensitive
titles ={'Diffusion','Imp/Exp','Diffusion + Imp/Exp'};
sens_sort = sort(sens,'descend');
sensi = find(sens>=sens_sort(15))';
for set = sensi
    xb = xb_orig(set,:);
    nn = 8;
    ps2ch = {[1 2 ], [4 7], [1 2 4 7]};
    del = linspace(0,-1,nn); % linspace(-3,3,nn);linspace(0,-2,nn); linspace(0,-2,nn); linspace(0,0.35,nn); linspace(0,-2,nn)];
    % del = repmat(-2:1:2,5,1);
    figure
    rgb = colormap('cool');
    rgb = rgb(round(linspace(20,length(rgb),nn))',:);
    
    n=1;
    

    
    for j = 1:length(ps2ch)
        xb_ = xb;
        for n = 1:nn
        xb_(ps2ch{j}) = xb(ps2ch{j}) + del(n);
        [e,p,protein,time,Beta,pnames,PROTEIN1X] = dsat_jp2(xb_);
        
        Y =  (linspace(0,1,51));
        Z =  (protein.nuclearDorsal.NC14(:,end));
        Zn = Z-min(Z(:)); Zn = Zn/max(Zn(:));
        subplot(2,3,j)
        hold on
        plot(Y,Z,'Color',rgb(n,:),'LineWidth',2);
        subplot(2,3,j+3)
        hold on
        plot(Y,Zn,'Color',rgb(n,:),'LineWidth',2);
        end
        subplot(2,3,j)
        title(titles{j})
    end
    %     view( [69 56])
    %     title(['\' pnames{i}],'FontSize',14)
    xlabel('value')
    %     xlabel('log_{10}(FC)')
    xlabel('DV Coordinate')
    ylabel('NC 14 gradient (Norm)')
%     delstr = num2cell(del);
%     for j = 1:nn
%         delstr{j} = num2str(delstr{j});
%     end
%     legend(delstr)
    
end