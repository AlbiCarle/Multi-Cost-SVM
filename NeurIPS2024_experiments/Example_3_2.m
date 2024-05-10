clc; clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data generation


mu1 = [4;6]; S1 = [1.3 0.9; 0.9 1.3];
mu2 = [3;8]; S2 = [.6 0; 0 1.4];

p_O = 0.0;

n = 10000;

X = []; Y = [];

p_A_array = linspace(0.1,0.9,9);

for p_A = p_A_array

    [X_p,Y_p] = generate_data_2(n/size(p_A_array,2), p_A, p_O, mu1, mu2, S1, S2);
    
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

%figure(1)

%% Exact Gaussian PSR

dimGrid=100; 

[K1, Z1] = meshgrid(linspace(min(Xts(:,1)), max(Xts(:,1)),dimGrid),...
                    linspace(min(Xts(:,2)), max(Xts(:,2)),dimGrid));

x=linspace(min(Xts(:,1)), max(Xts(:,1)), dimGrid);
y=linspace(min(Xts(:,2)), max(Xts(:,2)), dimGrid);
   
E=[K1(:) Z1(:)];

xmin = min(Xts(:,1));
xmax = max(Xts(:,1));

ymin = min(Xts(:,2));
ymax = max(Xts(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A1 = inv(S1);
A2 = inv(S2);
gamma = 0.5*log(det(S1)/det(S2)); 
Gamma_x = diag(0.5*(E'-mu1)'*A1*(E'-mu1)-0.5*(E'-mu2)'*A2*(E'-mu2)+gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
p_a1 = 0.2;  

rho = log((p_a1*epsilon)/((1-p_a1)*(1-epsilon)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target = -(Gamma_x-rho);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mycolor = [0,1,0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)

gscatter(Xts(:,1), Xts(:,2),Yts,'rb')

hold on

contourf(x, y, reshape(target,numel(y),numel(x)),[0.9999 0.9999] , ...
'linecolor', mycolor, 'LineWidth', 1,'FaceAlpha', 0.15);
legend('$\textrm{\textbf{x}}\hat{=}U$','$\textrm{\textbf{x}}\hat{=}S$','PSR surface', 'Interpreter','latex','Fontsize',16)
mycolormap = [mycolor(1)*ones(256,1), mycolor(2)*ones(256,1), mycolor(3)*ones(256,1)];

colormap(mycolormap)

xlim([xmin, xmax]);
ylim([ymin, ymax]);

y_pred_ts_gauss = -sign((diag(0.5*(Xts'-mu1)'*inv(S1)*(Xts'-mu1)-0.5*(Xts'-mu2)'*inv(S2)*(Xts'-mu2)+0.5*log(det(S1)/det(S2))))-rho);

[TPR_gau, FPR_gau, TNR_gau, FNR_gau, F1_gau, ACC_gau] = ConfusionMatrix(Yts, y_pred_ts_gauss,'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Training MCSVM

Tau  = rand(1,9); 
m = size(Tau,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kernel = 'polynomial';
param = 3;
eta = .001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha_bar = MCSVM_Train(Xtr, Ytr, kernel, param, Tau, eta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specialized on a dataset

n_ =  10000;
p_s = 0.2;

[X_p,Y_p] = generate_data_2(n_, p_s, p_O, mu1, mu2, S1, S2);
    
delta = 1E-4; 
epsilon = 0.1;

n_cl = ceil((7.47)/epsilon*log(1/delta));
r = ceil(epsilon*n_cl*0.5);
n_tr = 500;
n_ts = n_ - n_tr - n_cl -1;

[Xtr_p, Ytr_p, Xts_p, Yts_p, Xcl_p, Ycl_p] = split_dataset(X_p, Y_p, n_tr, n_ts, n_cl);

tau = 1-epsilon;

alpha_c = MCSVM_Train_c(Xtr, Ytr, Xcl_p, Ycl_p, kernel, param, tau, eta, alpha_bar);

b = offset_c(Xtr, Ytr, Xcl_p, Ycl_p, alpha_c, kernel, param, eta, tau, alpha_bar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)

gscatter(Xts_p(:,1), Xts_p(:,2),Yts_p,'rb')

hold on
plot_MCSVM(Xtr, Ytr, Xts_p, Yts_p, alpha_bar, b, 0, kernel, param, eta, mycolor,'--', 0.15);

xmin = min(Xts_p(:,1));
xmax = max(Xts_p(:,1));

ymin = min(Xts_p(:,2));
ymax = max(Xts_p(:,2));

xlim([xmin, xmax]);
ylim([ymin, ymax]);

legend('$\textrm{\textbf{x}}\hat{=}U$','$\textrm{\textbf{x}}\hat{=}S$','','Approximate PSR (cubic)', 'Interpreter','latex','Fontsize',16)

y_pred_ts = MCSVM_Test(Xtr, Ytr, Xts_p, alpha_bar, b, 0, kernel, param, eta);

[TPR_MCSVM, FPR_MCSVM, TNR_MCSVM, FNR_MCSVM, F1_MCSVM, ACC_MCSVM] = ConfusionMatrix(Yts_p, y_pred_ts,'on');

disp(['False positive rate:',num2str(FPR_MCSVM)])
