function [alpha] = MCSVM_Train_c(Xtr, Ytr, Xcl_p, Ycl_p, kernel, param, tau, eta, alpha_bar)

n = size(Xcl_p,1);

K = KernelMatrix(Xtr, Xcl_p, kernel, param);
D = diag(Ytr);
D_tilde = diag(Ycl_p);

P = alpha_bar'*(D*K*D_tilde);
f = -eta*P' + ones(n,1);
lb = zeros(n,1);
ub = 0.5*((1-2*tau).*Ycl_p + 1);
            
Aeq =  Ycl_p'; 
beq = 0;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXXXXXXXXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = linprog(-f,[],[],Aeq,beq,lb,ub);

%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXXXXXXXXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = x;


