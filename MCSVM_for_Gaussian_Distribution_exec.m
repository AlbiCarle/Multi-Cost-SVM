clc; clear all; close all;

% Configure it properly on your own device

addpath ./Utils/Algorithm/
addpath ./Utils/Evaluation_Visualization/
addpath ./Utils/Gaussian_Data_Generation/
addpath ./Utils/Various/

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

n_cl = 1000;
r = ceil(epsilon*n_cl*0.5);
n_tr = 3000;
n_ts = n - n_tr - n_cl -1;

[Xtr, Ytr, Xts, Yts, Xcl, Ycl] = split_dataset(X, Y, n_tr, n_ts, n_cl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)

gscatter(Xtr(:,1), Xtr(:,2), Ytr,'rb')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Exact Gaussian PSR

dimGrid=100; 

x=linspace(min(Xts(:,1)), max(Xts(:,1)), dimGrid);
y=linspace(min(Xts(:,2)), max(Xts(:,2)), dimGrid);

[K1, Z1] = meshgrid(x,y);
   
E=[K1(:) Z1(:)];

A1 = inv(S1); A2 = inv(S2);

gamma = 0.5*log(det(S1)/det(S2)); 
Gamma_x = diag(0.5*(E'-mu1)'*A1*(E'-mu1)-0.5*(E'-mu2)'*A2*(E'-mu2)+gamma);
   
p_S = 0.5;  

rho = log((p_S*epsilon)/((1-p_S)*(1-epsilon)));

target = -(Gamma_x-rho);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Training MCSVM

Tau  = rand(1,9); 
m = size(Tau,2);

kernel = 'polynomial';
param = 3;
eta = .001;

alpha_bar = MCSVM_Train(Xtr, Ytr, kernel, param, Tau, eta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specialized on a dataset

n_ =  10000;
p_s = 0.2;

[X_p,Y_p] = generate_data_2(n_, p_s, p_O, mu1, mu2, S1, S2);
    
n_cl = 1000;
r = ceil(epsilon*n_cl*0.5);
n_tr = 3000;
n_ts = n_ - n_tr - n_cl -1;

[Xtr_p, Ytr_p, Xts_p, Yts_p, Xcl_p, Ycl_p] = split_dataset(X_p, Y_p, n_tr, n_ts, n_cl);

%class ratio

tau = 1-epsilon;

alpha_c = MCSVM_Train_c(Xtr, Ytr, Xcl_p, Ycl_p, kernel, param, tau, eta, alpha_bar);

b = offset_c(Xtr, Ytr, Xcl_p, Ycl_p, alpha_c, kernel, param, eta, tau, alpha_bar);

mycolor = [0,1,0];

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

legend('$\textrm{\textbf{x}}\hat{=}U$','$\textrm{\textbf{x}}\hat{=}S$','','Approximate PSR (quintic)', 'Interpreter','latex','Fontsize',16)

saveas(gcf, 'Approx_cubic.fig');

cool_image_to_pdf('Approx_cubic', 'Approx_cubic')

y_pred_ts = MCSVM_Test(Xtr, Ytr, Xts_p, alpha_bar, b, 0, kernel, param, eta);

[TPR_MCSVM, FPR_MCSVM, TNR_MCSVM, FNR_MCSVM, F1_MCSVM, ACC_MCSVM] = ConfusionMatrix(Yts_p, y_pred_ts,'on');

disp(['False positive rate:',num2str(FPR_MCSVM)])
