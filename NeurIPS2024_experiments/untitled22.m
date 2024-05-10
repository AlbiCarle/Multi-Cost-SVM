clc; clear all; close all;
        
addpath ./SVM/

mu1 = [+1;+1]; S1m = [1.5 0.8;0.8 1.5];
mu2 = [-1;-1]; S2m = [2.5 0.1; 0.1 0.5];

p_S1 = 0.05;

p_O1 = 0.01;

n = 10000;

[X,Y] = generate_data_Albi(n, p_S1, p_O1, mu1, mu2, S1m, S2m);
    
[Z,C1,S1] = normalize(X,"norm",Inf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta = 1E-6; 
epsilon = 0.05;

n_cl = ceil((7.47)/epsilon*log(1/delta));
r = ceil(epsilon*n_cl*0.5);
n_tr = 1000;
n_ts = n - n_tr - n_cl -1;

[Xtr1, Ytr1, Xts1, Yts1, Xcl1, Ycl1] = split_dataset(Z, Y, n_tr, n_ts, n_cl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_S2 = 0.95;

p_O2 = 0.01;

[X,Y] = generate_data_Albi(n, p_S2, p_O2, mu1, mu2, S1m, S2m);
    
[Z,C2,S2] = normalize(X,"norm",Inf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta = 1E-6; 
epsilon = 0.05;

n_cl = ceil((7.47)/epsilon*log(1/delta));
r = ceil(epsilon*n_cl*0.5);
n_tr = 1000;
n_ts = n - n_tr - n_cl -1;

[Xtr2, Ytr2, Xts2, Yts2, Xcl2, Ycl2] = split_dataset(Z, Y, n_tr, n_ts, n_cl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kernel = 'linear';
param = .1;
eta = .01;

%%

tau = 0.5;

alpha1 = SSVM_Train(Xtr1.*S1+C1, Ytr1, kernel, param, tau, eta);
alpha2 = SSVM_Train(Xtr2.*S2+C2, Ytr2, kernel, param, tau, eta);

color1 = [255, 165, 0] / 255;  % Orange
color2 = [0, 128, 128] / 255;  % Teal

Xtrtilde1 = Xtr1.*S1+C1;
Xtstilde1 = Xts1.*S1+C1;

Xtrtilde2 = Xtr2.*S2+C2;
Xtstilde2 = Xts2.*S2+C2;

clf
figure('Position', [100, 100, 1000, 300]);

subplot(1, 3, 1);
% Plot SSVM
gscatter(Xtstilde1(:,1), Xtstilde1(:,2), Yts1, 'rb')
hold on
plot_SSVM2(Xtrtilde1, Ytr1, Xts1.*S1+C1, Yts1, alpha1, 0, 0, kernel, param, eta, color1,'-',0.5);


% Turn off legend
legend('off');

subplot(1, 3, 2);
gscatter(Xtstilde2(:,1), Xtstilde2(:,2), Yts2, 'rb')
hold on

plot_SSVM2(Xtr2.*S2+C2, Ytr2, Xts2.*S2+C2, Yts2, alpha2, 0, 0, kernel, param, eta, color2,'--',0.5);
% Turn off legend
legend('off');

% Set title
subplot(1, 3, 3);

plot_SSVM2(Xtr1.*S1+C1, Ytr1, Xts1.*S1+C1, Yts1, alpha1, 0, 0, kernel, param, eta, color1 ,'-',0);
hold on
plot_SSVM2(Xtr2.*S2+C2, Ytr2, Xts2.*S2+C2, Yts2, alpha2, 0, 0, kernel, param, eta, color2,'--',0);

tau2 = [0.4,0.7];

alpha12 = SSVM_Train(Xtr1, Ytr1, kernel, param, tau2, eta);
alpha22 = SSVM_Train(Xtr2, Ytr2, kernel, param, tau2, eta);

clf
figure('Position', [100, 100, 1000, 300]);

subplot(1, 3, 1);

% Plot SSVM
gscatter(Xtstilde1(:,1), Xtstilde1(:,2), Yts1, 'rb')
hold on
plot_SSVM2(Xtrtilde1, Ytr1, Xts1.*S1+C1, Yts1, alpha12, 0, 0, kernel, param, eta, color1,'-',0.5);

% Turn off legend
legend('off');

subplot(1, 3, 2);
gscatter(Xtstilde2(:,1), Xtstilde2(:,2), Yts2, 'rb')
hold on

plot_SSVM2(Xtr2.*S2+C2, Ytr2, Xts2.*S2+C2, Yts2, alpha22, 0, 0, kernel, param, eta, color2,'--',0.5);
% Turn off legend
legend('off');

% Set title
subplot(1, 3, 3);

plot_SSVM2(Xtr1.*S1+C1, Ytr1, Xts1.*S1+C1, Yts1, alpha12, 0, 0, kernel, param, eta, color1 ,'-',0);
hold on
plot_SSVM2(Xtr2.*S2+C2, Ytr2, Xts2.*S2+C2, Yts2, alpha22, 0, 0, kernel, param, eta, color2,'--',0);

tau5 = linspace(0.1,0.9,5);

alpha15 = SSVM_Train(Xtr1, Ytr1, kernel, param, tau5, eta);
alpha25 = SSVM_Train(Xtr2, Ytr2, kernel, param, tau5, eta);

clf
figure('Position', [100, 100, 1000, 300]);

subplot(1, 3, 1);
% Plot SSVM
gscatter(Xtstilde1(:,1), Xtstilde1(:,2), Yts1, 'rb')
hold on
plot_SSVM2(Xtrtilde1, Ytr1, Xts1.*S1+C1, Yts1, alpha15, 0, 0, kernel, param, eta, color1,'-',0.5);


% Turn off legend
legend('off');

subplot(1, 3, 2);
gscatter(Xtstilde2(:,1), Xtstilde2(:,2), Yts2, 'rb')
hold on

plot_SSVM2(Xtr2.*S2+C2, Ytr2, Xts2.*S2+C2, Yts2, alpha25, 0, 0, kernel, param, eta, color2,'--',0.5);
% Turn off legend
legend('off');

% Set title
subplot(1, 3, 3);

plot_SSVM2(Xtr1.*S1+C1, Ytr1, Xts1.*S1+C1, Yts1, alpha15, 0, 0, kernel, param, eta, color1 ,'-',0);
hold on
plot_SSVM2(Xtr2.*S2+C2, Ytr2, Xts2.*S2+C2, Yts2, alpha25, 0, 0, kernel, param, eta, color2,'--',0);

tau10 = linspace(0.1,0.9,10);

alpha110 = SSVM_Train(Xtr1, Ytr1, kernel, param, tau10, eta);
alpha210 = SSVM_Train(Xtr2, Ytr2, kernel, param, tau10, eta);

clf
figure('Position', [100, 100, 1000, 300]);

subplot(1, 3, 1);
% Plot SSVM
gscatter(Xtstilde1(:,1), Xtstilde1(:,2), Yts1, 'rb')
hold on
plot_SSVM2(Xtrtilde1, Ytr1, Xts1.*S1+C1, Yts1, alpha110, 0, 0, kernel, param, eta, color1,'-',0.2);

saveas(gcf, 'subplot110.fig'); % Save the first subplot


% Turn off legend
legend('off');
clf
figure('Position', [100, 100, 1000, 300]);

subplot(1, 3, 2);
gscatter(Xtstilde2(:,1), Xtstilde2(:,2), Yts2, 'rb')
hold on

plot_SSVM2(Xtr2.*S2+C2, Ytr2, Xts2.*S2+C2, Yts2, alpha210, 0, 0, kernel, param, eta, color2,'--',0.2);
% Turn off legend
legend('off');
saveas(gcf, 'subplot210.fig'); % Save the second subplot

clf
figure('Position', [100, 100, 1000, 300]);

% Set title
subplot(1, 3, 3);

plot_SSVM2(Xtr1.*S1+C1, Ytr1, Xts1.*S1+C1, Yts1, alpha110, 0, 0, kernel, param, eta, color1 ,'-',0);
hold on
plot_SSVM2(Xtr2.*S2+C2, Ytr2, Xts2.*S2+C2, Yts2, alpha210, 0, 0, kernel, param, eta, color2,'--',0);
