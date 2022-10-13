clk1 = clock;
clk = [num2str(clk1(1)),'-',num2strDU(clk1(2),2),'-',num2strDU(clk1(3),2),'_',...
	num2strDU(clk1(4),2),'-',num2strDU(clk1(5),2),'-',num2strDU(round(clk1(6)),2)];
id = strcat(clk,'_dsat_jp2_noshuttling0');

fhandle = 'dsat_jp2'; mm = 'min';
% lambda = 200; mu = max(ceil(lambda/7),41);

% lambda = 500; mu = 71;
% pf = .3; varphi = 1; G = 500; tmax=72*3600; % hrs * seconds/hr
lambda = 350; mu = 71;
pf = .3; varphi = 1; G = 50; tmax=72*3600; % hrs * seconds/hr

%      1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
% lu = [-1, 2,-4,-2,-5,-5,-3,-3,-2,-5,-5,-1,-3,-1,-5,-5, 2,-4,-5,-5;
%        0, 5, 5, 4,-1, 4, 4,-1, 5, 3, 4, 5, 3,-1,-1, 2, 5, 4,-1,-1];
%lu = [-1, 2,-5,-2,-5,-5,-4,-3,-3,-5,-5,-2,-2,-1,-5,-5, 1,-5,-5,-5;
%       1, 5, 5, 5,-1, 5, 5, 0, 5, 3, 4, 5, 4, 0,-1, 3, 5, 4,-1,-1];
% These limits are adapted from the previous successful runs

lu = [-4*ones(1,20); 4*ones(1,20)];
   
   
   
% lu = [-4*ones(1,20);2*ones(1,20)];
% lu(1,14)=-1;
% lu = [-3*ones(1,20); 2*ones(1,20)]; 
lu(1,14)=-1; lu(2,14)=0; % lower and upper limits
lu(1,2) = -12; lu(2,2) = -11;

[xb,Statistics,Gm]...
    =isres(fhandle,mm,lu,lambda,G,mu,pf,varphi,tmax);
% [xb,Statistics,F,X,Gm]...
%     =isres_mdo(fhandle,mm,lu,lambda,G,mu,pf,varphi,tmax);
save(id)

