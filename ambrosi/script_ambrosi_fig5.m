%script_ambrosi_fig5.m

load gyn.mat

COLR = [[255 211 201]/255;[212 175 167]/255;[158 104 92]/255;[99 48 37]/255];
% COLR = [[0 0 1];[1 0 0];[0 0.5 0];[1 0.5 0]];
% COLR = [[1 0 0];[1 0.5 0];[0.5 0.5 0];[0.2 0.5 0.2]];

[ gynr,~,~,~ ] = ambrosifun( El,Er,Eh,ntotal,R,S,xi,Gamma,ki,ke,Pcact,...
    kDeg,kb,r,n,tspan,Dl0,DlCact0,Cact0);
gyn = gynr / Dl0;
gynn = gyn(22:51,end); gynn = gynn-min(gynn(:)); %gynn = gynn/max(gynn(:));
gyn = gyn(:,end); gyn = gyn-min(gyn(:)); gyn = gyn/max(gyn(:));

Gamma = 2; 
[ gyn1r,~,~,~ ] = ambrosifun( El,Er,Eh,ntotal,R,S,xi,Gamma,ki,ke,Pcact,...
    kDeg,kb,r,n,tspan,Dl0,DlCact0,Cact0);
gyn1 = gyn1r / Dl0;
gyn1n = gyn1(22:51,end); gyn1n = gyn1n-min(gyn1n(:)); %gyn1n = gyn1n/max(gyn1n(:));
gyn1 = gyn1(:,end); gyn1 = gyn1-min(gyn1(:)); gyn1 = gyn1/max(gyn1(:));


ke = 1;
[ gyn2r,~,~,~ ] = ambrosifun( El,Er,Eh,ntotal,R,S,xi,Gamma,ki,ke,Pcact,...
    kDeg,kb,r,n,tspan,Dl0,DlCact0,Cact0);
gyn2 = gyn2r / Dl0;
gyn2n = gyn2(22:51,end); gyn2n = gyn2n-min(gyn2n(:)); %gyn2n = gyn2n/max(gyn2n(:));
gyn2 = gyn2(:,end); gyn2 = gyn2-min(gyn2(:)); gyn2 = gyn2/max(gyn2(:));


 Er = 117;
[ gyn3r,~,~,~ ] = ambrosifun( El,Er,Eh,ntotal,R,S,xi,Gamma,ki,ke,Pcact,...
    kDeg,kb,r,n,tspan,Dl0,DlCact0,Cact0);
gyn3 = gyn3r / Dl0;
gyn3n = gyn3(22:51,end); gyn3n = gyn3n-min(gyn3n(:)); %gyn3n = gyn3n/max(gyn3n(:));
gyn3 = gyn3(:,end); gyn3 = gyn3-min(gyn3(:)); gyn3 = gyn3/max(gyn3(:));

x = (1:length(gyn3))';

gynr = gynr(:,end);
gyn1r = gyn1r(:,end);
gyn2r = gyn2r(:,end);
gyn3r = gyn3r(:,end);

figure 
subplot(1,3,1)
hold on
plot(1:length(gynr),gynr,'. -','Color',COLR(1,:),'Markersize',12)
plot(1:length(gyn1r),gyn1r,'. -','Color',COLR(2,:),'Markersize',12)
plot(1:length(gyn2r),gyn2r,'. -','Color',COLR(3,:),'Markersize',12)
plot(1:length(gyn3r),gyn3r,'. -','Color',COLR(4,:),'Markersize',12)
% scatter(1:length(gynr),gynr,[],COLR(1,:),'Filled')
% scatter(1:length(gyn1r),gyn1r,[],COLR(2,:),'Filled')
% scatter(1:length(gyn2r),gyn2r,[],COLR(3,:),'Filled')
% scatter(1:length(gyn3r),gyn3r,[],COLR(4,:),'Filled')
title('Raw')

subplot(1,3,2)
hold on
plot(x,gyn,'. -','Color',COLR(1,:),'Markersize',12)
plot(x,gyn1,'. -','Color',COLR(2,:),'Markersize',12)
plot(x,gyn2,'. -','Color',COLR(3,:),'Markersize',12)
plot(x,gyn3,'. -','Color',COLR(4,:),'Markersize',12)
% scatter(x,gyn,[],COLR(1,:),'Filled')
% scatter(x,gyn1,[],COLR(2,:),'Filled')
% scatter(x,gyn2,[],COLR(3,:),'Filled')
% scatter(x,gyn3,[],COLR(4,:),'Filled')
title('normalized the right way')


subplot(1,3,3)
hold on
plot(1:30,gynn,'. -','Color',COLR(1,:),'Markersize',12)
beta = sum(gyn1n(:).*gynn(:))/sum(gyn1n(:).^2);
plot(1:30,beta*gyn1n,'. -','Color',COLR(2,:),'Markersize',12)
beta = sum(gyn2n(:).*gynn(:))/sum(gyn2n(:).^2);
plot(1:30,beta*gyn2n,'. -','Color',COLR(3,:),'Markersize',12)
beta = sum(gyn3n(:).*gynn(:))/sum(gyn3n(:).^2);
plot(1:30,beta*gyn3n,'. -','Color',COLR(4,:),'Markersize',12)
% scatter(1:30,gynn,[],COLR(1,:),'Filled')
% beta = sum(gyn1n(:).*gynn(:))/sum(gyn1n(:).^2);
% scatter(1:30,beta*gyn1n,[],COLR(2,:),'Filled')
% beta = sum(gyn2n(:).*gynn(:))/sum(gyn2n(:).^2);
% scatter(1:30,beta*gyn2n,[],COLR(3,:),'Filled')
% beta = sum(gyn3n(:).*gynn(:))/sum(gyn3n(:).^2);
% scatter(1:30,beta*gyn3n,[],COLR(4,:),'Filled')

title('normalized the other way')