clc; clear all; close all;

% Data generation

%addpath /Users/albertocarlevaro/Documents/Albi/Fabrizio&Teo/'Probabilistic Safety regions'/'Conformal Scalable Classifiers'/SVM
addpath ./SVM
addpath /Users/albertocarlevaro/Documents/Albi/Fabrizio&Teo/'Probabilistic Safety regions'/'SCforExpFamilies'/Exponential_Distributions/

mu1 = [4;6]; S1 = [1.3 0.9; 0.9 1.3];
mu2 = [3;8]; S2 = [.6 0; 0 1.4];

p_O = 0.0;

n = 10000;

X = []; Y = [];

p_A_array = linspace(0.1,0.9,9);

for p_A = p_A_array

    [X_p,Y_p] = generate_data_Albi(n/size(p_A_array,2), p_A, p_O, mu1, mu2, S1, S2);
    
    X = [X; X_p]; Y = [Y;Y_p];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta = 1E-3; 
epsilon = 0.05;

n_cl = ceil((7.47)/epsilon*log(1/delta));
r = ceil(epsilon*n_cl*0.5);
n_tr = 500;
n_ts = n - n_tr - n_cl -1;

[Xtr, Ytr, Xts, Yts, Xcl, Ycl] = split_dataset(X, Y, n_tr, n_ts, n_cl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)

gscatter(Xtr(:,1), Xtr(:,2), Ytr,'rb')

%% Exact Gaussian PSR

dimGrid=100; 

[K1, Z1] = meshgrid(linspace(min(Xts(:,1)), max(Xts(:,1)),dimGrid),...
                    linspace(min(Xts(:,2)), max(Xts(:,2)),dimGrid));

x=linspace(min(Xts(:,1)), max(Xts(:,1)), dimGrid);
y=linspace(min(Xts(:,2)), max(Xts(:,2)), dimGrid);
   
E=[K1(:) Z1(:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A1 = inv(S1);
A2 = inv(S2);
gamma = 0.5*log(det(S1)/det(S2)); 
Gamma_x = diag(0.5*(E'-mu1)'*A1*(E'-mu1)-0.5*(E'-mu2)'*A2*(E'-mu2)+gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
p_a1 = 0.5;  

rho = log((p_a1*epsilon)/((1-p_a1)*(1-epsilon)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target = -(Gamma_x-rho);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mycolor = [0,1,1];

figure(2)

gscatter(Xts(:,1), Xts(:,2),Yts,'rb')

hold on

contourf(x, y, reshape(target,numel(y),numel(x)),[0.9999 0.9999] , ...
'linecolor', mycolor, 'LineWidth', 1,'FaceAlpha', 0.15);
legend({'$-1$','$+1$','Gaussian PSR'},'Interpreter','latex')
mycolormap = [mycolor(1)*ones(256,1), mycolor(2)*ones(256,1), mycolor(3)*ones(256,1)];

colormap(mycolormap)

y_pred_ts_gauss = -sign((diag(0.5*(Xts'-mu1)'*inv(S1)*(Xts'-mu1)-0.5*(Xts'-mu2)'*inv(S2)*(Xts'-mu2)+0.5*log(det(S1)/det(S2))))-rho);


[TPR_gau, FPR_gau, TNR_gau, FNR_gau, F1_gau, ACC_gau] = ConfusionMatrix(Yts, y_pred_ts_gauss,'on');
%% Training SSVM

%tau = linspace(0.1,0.9,10);
Tau  = [0.1,0.5,0.9];%rand(1,2); 
m = size(Tau,2);

kernel = 'linear';
param = 2;
eta = .01;

[alpha_bar, alpha_k] = SSVM_Train(Xtr, Ytr, kernel, param, Tau, eta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = {};

for k = 1:m
    C{k} = offset(Xtr, Ytr, alpha_k{k}, kernel, param, eta, Tau(k));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_pred_ts = SSVM_Test(Xtr, Ytr, Xts, alpha_bar, 0, 0, kernel, param, eta);

[TPR_SSVM, FPR_SSVM, TNR_SSVM, FNR_SSVM, F1_SSVM, ACC_SSVM] = ConfusionMatrix(Yts, y_pred_ts,'on');

figure(3)

gscatter(Xts(:,1), Xts(:,2),Yts,'rb')

hold on
plot_SSVM2(Xtr, Ytr, Xts, Yts, alpha_bar, 0, 0, kernel, param, eta, mycolor,'--', 0.01);
legend({'$-1$','$+1$',' ', 'Approximate PSR'},'Interpreter','latex')

%% Hold on the two regions

figure(4)

gscatter(Xts(:,1), Xts(:,2),Yts,'rb')

hold on
plot_SSVM2(Xtr, Ytr, Xts, Yts, alpha_bar, 0, 0, kernel, param, eta, [0,1,0],'--', 0.15);
legend({'$-1$','$+1$',' ', 'Approximate PSR'},'Interpreter','latex')

hold on

contourf(x, y, reshape(target,numel(y),numel(x)),[0.9999 0.9999] , ...
'linecolor', mycolor, 'LineWidth', 2,'FaceAlpha', 0.15);
legend({'$-1$','$+1$','Gaussian PSR'},'Interpreter','latex')
mycolormap = [mycolor(1)*ones(256,1), mycolor(2)*ones(256,1), mycolor(3)*ones(256,1)];

colormap(mycolormap)

%% Calibrate the SSVM with scaling

Xcl_U = Xcl(Ycl == -1,:); Ycl_U = Ycl(Ycl == -1);
Xcl_S = Xcl(Ycl == +1,:); Ycl_S = Ycl(Ycl == +1);

n_U = size(Xcl_U,1);
n_S = size(Xcl_S,1);

n_c = n_S + n_U ;

r = ceil(epsilon*n_c*0.5);

Gamma_rho = barrhoSVM(Xtr, Ytr, Xcl_U, 0, alpha_bar, kernel, param, eta);

Gamma_rho_sort = sort(Gamma_rho,'descend');

rho_star = Gamma_rho_sort(r);

y_pred_ts_sc = SSVM_Test(Xtr, Ytr, Xts, alpha_bar, 0, rho_star, kernel, param, eta);

[TPR_SSVM_sc, FPR_SSVM_sc, TNR_SSVM_sc, FNR_SSVM_sc, F1_SSVM_sc, ACC_SSVM_sc] = ConfusionMatrix(Yts, y_pred_ts_sc,'on');

figure(5)

gscatter(Xts(:,1), Xts(:,2),Yts,'rb')

hold on
plot_SSVM2(Xtr, Ytr, Xts, Yts, alpha_bar, 0, rho_star, kernel, param, eta, [0,0,0],'--', 0.1);

hold on
plot_SSVM2(Xtr, Ytr, Xts, Yts, alpha_bar, 0, 0, kernel, param, eta, [0,1,0],'--', 0.1);

hold on

contourf(x, y, reshape(target,numel(y),numel(x)),[0.9999 0.9999] , ...
'linecolor', mycolor, 'LineWidth', 2,'FaceAlpha', 0.1);
legend({'$-1$','$+1$','Gaussian PSR'},'Interpreter','latex')
mycolormap = [mycolor(1)*ones(256,1), mycolor(2)*ones(256,1), mycolor(3)*ones(256,1)];

colormap(mycolormap)

legend({'$-1$','$+1$','','Scaling SSVM','', 'SSVM', 'True PSR'},'Interpreter','latex')

%%
disp(['$\varepsilon$ = ',num2str(epsilon)])
disp(['FPR Gaussian PSR = ', num2str(FPR_gau)]);
disp(['FPR SSVM PSR = ', num2str(FPR_SSVM)]);
disp(['FPR SSVM with scaling PSR = ', num2str(FPR_SSVM_sc)]);

