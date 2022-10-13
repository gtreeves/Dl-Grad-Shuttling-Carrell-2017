function [ un,uc,wc,vc ] = ambrosifun( El,Er,Eh,ntotal,R,S,xi,Gamma,ki,ke,Pcact,kDeg,kb,r,n,tspan,Dl0,DlCact0,Cact0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Vcl = 4*pi/3*El*Er*Er-((4*pi)/3)*(El - Eh)*(Er - Eh)*(Er - Eh) ;
Vcomp = Vcl/ntotal;
L = 2*pi*Er;
T = n*Vcomp/(pi*Er*Er-pi*(Er-Eh)^2);
Am = Eh*T;
Vn = 4/3*pi*r^3;
An = 4*pi*r^2;
Vc = Vcomp - Vn;
h = 1:n;
x = L*((h-0.5)/n - 0.5);
Kd = R./(S+abs(x).^xi)';

e = ones(n,1); P = spdiags([e -2*e e],[-1 0 1],n,n);
P(1,end)=1; P(end,1) = 1;

[t,y] = ode15s(@ambfun,[0 tspan],[zeros(n,1);Dl0*ones(n,1);DlCact0*ones(n,1);Cact0*ones(n,1)] ,[],n,An,ki,Vn,ke,Gamma,Am,Vc,P,Kd,kb,Pcact,kDeg);
    y=y';
    un = y(1:n,:);
    uc = y(n+1:2*n,:);
    wc = y(2*n+1:3*n,:);
    vc = y(3*n+1:4*n,:);
end

function dydt = ambfun(~,y,n,An,ki,Vn,ke,Gamma,Am,Vc,P,Kd,kb,Pcact,kDeg)
    un = y(1:n);
    uc = y(n+1:2*n);
    wc = y(2*n+1:3*n);
    vc = y(3*n+1:4*n);
    
    if length(Gamma)==1
        Gamma = Gamma*ones(3,1);
    end
    
    f1 = An*ki/Vn*uc-An*ke/Vn*un;
    f2 = Gamma(1)*Am/Vc*P*uc+Kd.*wc-kb*uc.*vc-(An*ki/Vc*uc-An*ke/Vc*un);
    f3 = Gamma(2)*Am/Vc*P*wc-Kd.*wc+kb*uc.*vc;
    f4 = Gamma(3)*Am/Vc*P*vc+Kd.*wc-kb*uc.*vc+Pcact-kDeg*vc;
    
    dydt = [f1;f2;f3;f4];
end