clc; clear all; close all;

% Configure it properly on your own device

addpath ./Utils/Algorithm/
addpath ./Utils/Evaluation_Visualization/
addpath ./Utils/Gaussian_Data_Generation/
addpath ./Utils/Various/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu1 = [4;6]; S1 = [1.3 0.9; 0.9 1.3];
mu2 = [3;8]; S2 = [.6 0; 0 1.4];

p_O = 0.0;

n = 5000;

X = []; Y = [];

p_A_array = [0.1, 0.2, 0.5, 0.9];

epsilon = 0.05;

mycolor =  [0, 1, 0];

i = 0;

figure('Position', [100, 100, 1000, 800]); % Adjust the position and size as needed

for p_a = p_A_array

    i = i+1;

    [X_p,Y_p] = generate_data_2(n, p_a, p_O, mu1, mu2, S1, S2);
    
    X = [X; X_p]; Y = [Y;Y_p];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    epsilon = 0.05;
    
    n_cl = 1;
    n_tr = 1;
    n_ts = n - n_tr - n_cl -1;
    
    [Xtr, Ytr, Xts, Yts, Xcl, Ycl] = split_dataset(X, Y, n_tr, n_ts, n_cl);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xmin = min(Xts(:,1));
    xmax = max(Xts(:,1));
    
    ymin = min(Xts(:,2));
    ymax = max(Xts(:,2));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dimGrid=50; 

    x=linspace(min(Xts(:,1)), max(Xts(:,1)), dimGrid);
    y=linspace(min(Xts(:,2)), max(Xts(:,2)), dimGrid);
    
    [K1, Z1] = meshgrid(linspace(min(Xts(:,1))-1, max(Xts(:,1))+1,dimGrid),...
                        linspace(min(Xts(:,2))-1, max(Xts(:,2))+1,dimGrid));
    
    E=[K1(:) Z1(:)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A1 = inv(S1);
    A2 = inv(S2);

    gamma = 0.5*log(det(S1)/det(S2)); 
    Gamma_x = diag(0.5*(E'-mu1)'*A1*(E'-mu1)-0.5*(E'-mu2)'*A2*(E'-mu2)+gamma);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rho = log((p_a*epsilon)/((1-p_a)*(1-epsilon)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    target = -(Gamma_x-rho);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    subplot(4, 4, i);

    gscatter(Xts(:,1), Xts(:,2),Yts,'rb')
   
    hold on
    
    contourf(x, y, reshape(target,numel(y),numel(x)),[0.9999 0.9999] , ...
    'linecolor', mycolor, 'LineWidth', 1,'FaceAlpha', 0.15);
    legend({'$-1$','$+1$','Gaussian PSR'},'Interpreter','latex')
    mycolormap = [mycolor(1)*ones(256,1), mycolor(2)*ones(256,1), mycolor(3)*ones(256,1)];

    colormap(gca, mycolormap);
    
    xlim([xmin, xmax]);
    ylim([ymin, ymax]);
    legend off
    hold off
    title(['$p_S \; = \; $', num2str(p_a)],'Interpreter','latex','FontSize',16);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Epsilon = [0.01, 0.05, 0.1, 0.5];

P_A = 0.5;

[X_p,Y_p] = generate_data_2(n, P_A, p_O, mu1, mu2, S1, S2);
    
X = [X; X_p]; Y = [Y;Y_p];

j = i;

for epsilon = Epsilon

    j = j+1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    n_cl = 1;
    n_tr = 1;
    n_ts = n - n_tr - n_cl -1;
    
    [Xtr, Ytr, Xts, Yts, Xcl, Ycl] = split_dataset(X, Y, n_tr, n_ts, n_cl);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xmin = min(Xts(:,1));
    xmax = max(Xts(:,1));
    
    ymin = min(Xts(:,2));
    ymax = max(Xts(:,2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dimGrid=50; 
    
    [K1, Z1] = meshgrid(linspace(min(Xts(:,1))-1, max(Xts(:,1))+1,dimGrid),...
                        linspace(min(Xts(:,2))-1, max(Xts(:,2))+1,dimGrid));
    
    x=linspace(min(Xts(:,1)), max(Xts(:,1)), dimGrid);
    y=linspace(min(Xts(:,2)), max(Xts(:,2)), dimGrid);
       
    E=[K1(:) Z1(:)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    A1 = inv(S1);
    A2 = inv(S2);
    
    gamma = 0.5*log(det(S1)/det(S2)); 
    Gamma_x = diag(0.5*(E'-mu1)'*A1*(E'-mu1)-0.5*(E'-mu2)'*A2*(E'-mu2)+gamma);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rho = log((P_A*epsilon)/((1-P_A)*(1-epsilon)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    target = -(Gamma_x-rho);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    subplot(4, 4, j);

    gscatter(Xts(:,1), Xts(:,2),Yts,'rb')
   
    hold on
    
    contourf(x, y, reshape(target,numel(y),numel(x)),[0.9999 0.9999] , ...
    'linecolor', mycolor, 'LineWidth', 1,'FaceAlpha', 0.15);
    legend({'$-1$','$+1$','Gaussian PSR'},'Interpreter','latex')
    mycolormap = [mycolor(1)*ones(256,1), mycolor(2)*ones(256,1), mycolor(3)*ones(256,1)];

    colormap(gca, mycolormap);
    
    xlim([xmin, xmax]);
    ylim([ymin, ymax]);
    legend off
    hold off
    title(['$\varepsilon \; = \; $', num2str(epsilon)],'Interpreter','latex','FontSize',16);

end

subplot(4,4,[9 10 11 12 13 14 15 16])

dimGrid=25; 

[K1, Z1] = meshgrid(linspace(min(Xts(:,1))-1, max(Xts(:,1))+1,dimGrid),...
                    linspace(min(Xts(:,2))-1, max(Xts(:,2))+1,dimGrid));

x=linspace(min(Xts(:,1)), max(Xts(:,1)), dimGrid);
y=linspace(min(Xts(:,2)), max(Xts(:,2)), dimGrid);
   
E=[K1(:) Z1(:)];

A1 = inv(S1);
A2 = inv(S2);

gamma = 0.5*log(det(S1)/det(S2)); 
Gamma_x = diag(0.5*(E'-mu1)'*A1*(E'-mu1)-0.5*(E'-mu2)'*A2*(E'-mu2)+gamma);

grid on

gscatter(Xts(:,1), Xts(:,2),Yts,'rb')
   
hold on

sc = meshc(K1,Z1, reshape(-Gamma_x, numel(x), numel(y)),'LineWidth',2,'FaceAlpha', 0.5);
ax = gca;
%ax.ZLim(2) = 70;
%sc(2).ZLocation = 'zmax';
colormap(gca, mycolor);
title('PSR decision surface ($\rho(p_S,\varepsilon) = 0$)','Interpreter','latex','FontSize',16);
legend('$\textrm{\textbf{x}}\hat{=}U$','$\textrm{\textbf{x}}\hat{=}S$','PSR surface', 'Interpreter','latex','Fontsize',16)
