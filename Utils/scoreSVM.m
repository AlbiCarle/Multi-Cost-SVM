function s = scoreSVM(Xtr, Ytr, X, Y, b, alpha, kernel, param, eta)

if isequal(kernel, 'linear')

     w = -eta*(alpha'*diag(Ytr))*Xtr;

     r = b - X*w';

else
    
    K = KernelMatrix(Xtr,X,kernel,param);
    D = diag(Ytr);
    
    r = b + eta*(K'*(D*alpha));
end

    s = -Y.*r;
