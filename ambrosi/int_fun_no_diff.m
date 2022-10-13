function F = int_fun_no_diff(~,y,M,sigmaU,muU,...
    gamma,psi,alpha,beta,phi,xi,Vn,Vc,An)
if M>1
    x=linspace_offset(0,1,M)';
else
    x = 0;
end
un = y(1:M);
uc = y(M+1:2*M);
wc = y(2*M+1:3*M);
vc = y(3*M+1:4*M);

f1 = (sigmaU*An*uc-muU*An*un)/Vn;
f2 = (Vc*(beta./(phi+x.^xi)).*wc-gamma*Vc*(uc.*vc)...
    -sigmaU*An*uc+muU*An*un)/Vc;
f3 = (-Vc*(beta./(phi+x.^xi)).*wc+gamma*Vc*(uc.*vc))/Vc;
f4 = (psi*Vc*(beta./(phi+x.^xi)).*wc-gamma*psi*Vc*(uc.*vc)+...
    1-alpha*Vc*vc)/Vc;

F = [f1;f2;f3;f4];