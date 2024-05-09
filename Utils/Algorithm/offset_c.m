function b = offset_c(Xtr, Ytr,Xcl_p, Ycl_p, alpha, kernel, param, eta, tau, alpha_bar)

C = 0.5*((1-2*tau).*Ycl_p + 1);
thr = 0.001;

ind = find(alpha-(0+thr)>0 & alpha-(C-thr)<0);

X_SV = Xcl_p(ind,:); Y_SV = Ycl_p(ind,:); 

if isequal(kernel, 'linear')

    w = -eta*(alpha_bar'*diag(Ytr))*Xtr;

    bii = X_SV*w' + Y_SV;

else

    K = KernelMatrix(Xtr,X_SV,kernel,param);
    D =  diag(Ytr);

    bii = -eta*(K'*(D*alpha_bar)) +Y_SV;
    
    
end

b = randsample(bii,1);
