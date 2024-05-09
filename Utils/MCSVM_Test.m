function y = MCVM_Test(Xtr, Ytr, Xts, alpha_bar, b, b_eps, kernel, param, eta)

K = KernelMatrix(Xtr,Xts,kernel,param);
D = diag(Ytr);

y = -sign(-b -eta*(K'*(D*alpha_bar)) + b_eps); 
