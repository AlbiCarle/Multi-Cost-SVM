function b = offset(Xtr, Ytr, alpha, kernel, param, eta, tau)

C = 0.5*((1-2*tau).*Ytr + 1);
thr = 0.00001;

ind = find(alpha-(0+thr)>0 & alpha-(C-thr)<0);

X_SV = Xtr(ind,:); Y_SV = Ytr(ind,:); 

if isequal(kernel, 'linear')

    w = -eta*(alpha'*diag(Ytr))*Xtr;

    bii = X_SV*w' + Y_SV;

else

    K = KernelMatrix(Xtr,X_SV,kernel,param);
    D =  diag(Ytr);

    bii = -eta*(K'*(D*alpha)) +Y_SV;
    
end

b = randsample(bii,1);
