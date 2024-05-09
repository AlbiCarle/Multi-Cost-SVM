function [] = plot_SSVM(Xtr, Ytr, Xts, Yts, alpha_bar, b, b_eps, kernel, param, eta, mycolor,linesty,gsc)

ax = gca;
        
xlim(ax.XLim); % Set the x-axis limits
ylim(ax.YLim); % Set the y-axis limits

x_limits = ax.XLim;
y_limits = ax.YLim;

% Calculate dimGrid based on the range of x and y
x_range = diff(x_limits);
y_range = diff(y_limits);

dimGrid=100; % dimGrid*dimGrid

% Adjust dimGrid based on the ranges
if x_range > y_range
    dimGrid = round(dimGrid * x_range / y_range);
else
    dimGrid = round(dimGrid * y_range / x_range);
end

[K1, Z1] = meshgrid(linspace(x_limits(1) - 1, x_limits(2) + 1, dimGrid), ...
                    linspace(y_limits(1) - 1, y_limits(2) + 1, dimGrid));

x = linspace(x_limits(1) - 1, x_limits(2) + 1, dimGrid);
y = linspace(y_limits(1) - 1, y_limits(2) + 1, dimGrid);
   
K1=K1(:); Z1=Z1(:);
E=[K1 Z1];

if isequal(kernel, 'linear')

    %close all
    w = -eta*(alpha_bar'*diag(Ytr))*Xtr;

    f = @(z) 1/w(2)*((b-b_eps) - w(1)*z);
    y_lin = f(x);

    % Plot the contour of w^Tx - b = 0 (Hyperplane)
    %contour(x, y, func_values_matrix, [0 0], 'LineColor', mycolor, 'LineWidth', 2);
    %hold on;
    hyperplane_func = @(t) w*t' -b +b_eps;
    func_values = hyperplane_func(E);
    func_values_matrix = reshape(func_values, length(x), length(y));
    filled_region = func_values_matrix <= 0;
    contourf(x, y, filled_region, [1 1], 'LineColor', 'none', 'FaceColor', mycolor, 'FaceAlpha', 0.2);
   
    hold on

    plot(x,y_lin,'Color',mycolor,'LineStyle',linesty,'LineWidth',2)

    legend off

    %set(gca,'xtick',[])
    %    set(gca,'xticklabel',[])
    %    set(gca,'ytick',[])
    %    set(gca,'yticklabel',[])
    %   axis off
    %axis equal
    
else

    y_pred = SSVM_Test(Xtr, Ytr, E, alpha_bar, b, b_eps, kernel, param, eta);

  
    contourf(x, y, reshape(y_pred,numel(y),numel(x)),[0.9999 0.9999] , ...
    'linecolor', mycolor, 'LineWidth', 2,'FaceAlpha', 0.15,'LineStyle',linesty);
    mycolormap = [mycolor(1)*ones(256,1), mycolor(2)*ones(256,1), mycolor(3)*ones(256,1)];

    colormap(mycolormap)
    hold off

end


      % %gscatter(Xvl(:,1), Xvl(:,2), Yvl, 'rb')
        % 
        % 
        % % Set the x-axis and y-axis limits
        % xlim([ax.XLim(1), ax.XLim(2)]); % Set the limits for the x-axis
        % ylim([ax.YLim(1), ax.YLim(2)]); % Set the limits for the y-axis
        % 
        %  corners = [ax.XLim(1), ax.YLim(1);...
        %             ax.XLim(2), ax.YLim(2);...
        %            ax.XLim(1), ax.YLim(2);...
        %            ax.XLim(2), ax.YLim(1)];
        % 
        % a = sign(w*corners');
        % corn_idxs = find(a>0);
        % 
        % m = -w(1)/w(2);
        % 
        % if m <0
        %     corner1 = corners(corn_idxs(1),:);
        %     corner2 = corners(corn_idxs(2),:);
        % elseif m > 0
        %     corner1 = corners(corn_idxs(2),:);
        %     corner2 = corners(corn_idxs(1),:);
        % end
        % 
        % hold on
        % 
        % fill([x(dimGrid) corner1(1) corner2(1) x(1)],[f(x(dimGrid)) corner2(2) corner1(2) f(x(1))], mycolor,'LineStyle','none',...
        %    'FaceAlpha',0.05, 'EdgeColor', mycolor, 'LineWidth',0.001);
        % 
        %  hold on
        % 
        % line([x(1), corner1], [f(x(1)), corner2], 'Color', [1,1,1], 'LineWidth', 1);
        % line([x(1), x(dimGrid)], [f(x(1)), f(x(dimGrid))], 'Color', [1,1,1], 'LineWidth', 1);
        % line([x(dimGrid), x(dimGrid)], [f(x(dimGrid)), corner2], 'Color', [1,1,1], 'LineWidth', 1);
        % 
        % hold on