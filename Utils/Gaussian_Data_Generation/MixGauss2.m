function [X, Y] = MixGauss2(mu, Sigma, n)

d = size(mu,1);
p = size(mu,2);

X = [];
Y = [];

for i = 1:p
    Xi = zeros(n,d);
    Yi = zeros(n,1);
    for j = 1:n
        x = mvnrnd(mu,Sigma,n); 
        Xi(j,:) = x;
        Yi(j) = i-1;
    end
    X = [X; Xi];
    Y = [Y; Yi];
end
