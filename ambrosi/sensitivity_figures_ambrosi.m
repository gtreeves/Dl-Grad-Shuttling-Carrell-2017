%% sensitivity analysis
load('results999999_dsat_jp2_mdo.mat','xb')
% close all
xb_orig = xb;
b = 1;
for set = 6:10
xb = xb_orig(set,:);
nn = 8;
ps2ch = [4 7 1 2 14 ];
del = [linspace(-3,3,nn); linspace(-3,3,nn);linspace(0,-2,nn); linspace(0,-2,nn); linspace(0,0.35,nn); linspace(0,-2,nn)];
% del = repmat(-2:1:2,5,1);
figure
rgb = colormap('cool');
rgb = rgb(round(linspace(20,length(rgb),nn))',:);

n=1;
for i = ps2ch
    xb_ = xb;
    subplot(3,2,n)
     hold on;
%     if i==14
%         xb_(1,13) = xb_(1,13)+1; %THIS INCREASES BETA WHEN MESSING WITH PHI
%     end
    for j = 1:nn
        xb_(1,i) = xb(1,i) + del(n,j);
        XB{b,1}=xb_; b=b+1;
        [e,p,protein,time,Beta,pnames,PROTEIN1X] = dsat_jp2(xb_);
        X = del(n,j)*ones(1,102);
%         X = xb_(1,i)*ones(1,102);
        Y =  (linspace(0,1,51));
        Z =  (protein.nuclearDorsal.NC14(:,end));
        Z = Z-min(Z(:)); Z = Z/max(Z(:));
        plot(Y,Z,'Color',rgb(j,:),'LineWidth',2);
    end
%     view( [69 56])
    title(['\' pnames{i}],'FontSize',14)
    xlabel('value')
%     xlabel('log_{10}(FC)')
    xlabel('DV Coordinate')
    ylabel('NC 14 gradient (Norm)')
    delstr = num2cell(del(n,:));
    for j = 1:nn
        delstr{j} = num2str(delstr{j});
    end
    legend(delstr)
    n=n+1;
end

    xb_ = xb;
    subplot(3,2,n)
     hold on;
    Ddose = linspace(0.5,1,nn);
    for j = 1:nn
        
        [e,p,protein,time,Beta,pnames,PROTEIN1X] = dsat_jp2(xb_,'dl',Ddose(j));
        X = 2*Ddose(j)*ones(1,102);
        Y = (linspace(0,1,51));
        Z =  (protein.nuclearDorsal.NC14(:,end));
        Z = Z-min(Z(:)); Z = Z/max(Z(:));
        plot(Y,Z,'Color',rgb(j,:),'LineWidth',2);
    end
%     view( [69 56])
    title('dl dosage','FontSize',14)
%     xlabel('value')
    xlabel('DV Coordinate')
    ylabel('NC 14 gradient (Norm)')
    n=n+1;
    
end
