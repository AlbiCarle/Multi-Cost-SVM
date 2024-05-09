function [alpha_bar, alpha_k] = SSVM_Train(Xtr, Ytr, kernel, param, tau, eta)

n = size(Xtr,1);
if(size(tau,1)>1)
    tau = tau';
end

m = size(tau,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXXXXXXXXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%tau = sort(tau,'descend');

%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXXXXXXXXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = KernelMatrix(Xtr, Xtr, kernel, param);
D = diag(Ytr);

%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXXXXXXXXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x = [alpha_bar_1, alpha_bar_2, ..., alpha_bar_n, -------------> n var
%      alpha_{1,1}, alpha_{1,2}, ..., alpha_{1,m}, --> m var}       +
%      alpha_{2,1}, alpha_{2,2}, ..., alpha_{2,m}, --> m var} --> nm var
%      ...    ...     ...     ...     ...     ...,          }       |
%      alpha_{n,1}, alpha_{n,2}, ..., alpha_{n,m}] --> m var}       |
%                                                                 n+nm = n(m+1) var                                                              

% H = [H1 0^T; 0 0^T] where
% H1 = D*K*D, size(H1) = [n,n]
% size(0) = [nm,n]

%%%%%                   %%%%
%%%%%   min x'Hx +f'x   %%%%
%%%%%   s.t.            %%%%

% min x^THx +\Sum{i=1}{n}alpha_bar_i
% s.t.
% {alpha_bar_i = \Sum{k=1}{m} alpha_{i,k}
%  --> alpha_bar_i - \Sum{k=1}{m} alpha_{i,k} = 0
%       --> alpha_bar_1 - \Sum{k=1}{m} alpha_{1,k} = 0
%           --> Aeq_1(1,:) = [1,0, ..., 0, 1, ..., 1, 0, ..., 0, ..., 0, ..., 0]
%       --> alpha_bar_2 - \Sum{k=1}{m} alpha_{2,k} = 0
%           --> Aeq_1(2,:) = [0,1, ..., 0, 0, ..., 0, 1, ..., 1, ..., 0, ..., 0]
%       ...     ...     ...     ...     ...     ...     
%       --> alpha_bar_n - \Sum{k=1}{m} alpha_{n,k} = 0
%           --> Aeq_1(n,:) = [0,0, ..., 1, 0, ..., 0, 0, ..., 0, ..., 1, ..., 1]
%
% {\Sum{i=1}{n} y_i*alpha_{i,k} = 0
% {0 <= \alpha_{i,k} <= 1/2*((1-2\tau_k)y_i +1)

% H
disp('ehy ehy')
            H1 = 1*eta*D*K*D;
            O = sparse(zeros(n,n*m));
            OO = sparse(zeros(n*m,n*m));
            
            H = [H1 O ; O' OO]; %remember, in quadprog matlab there is already 1/2
            
            %f = ones(n*(m+1),1); WARNING, YOU CHANGED f
            f = [ones(n,1);zeros(n*m,1)];
            
            % unequality constraints
            
            lb = zeros(n*(m+1),1);
            ub = ones(n*(m+1),1);
            
            %tau_bar = sum(tau);
            
            %C = m/2*((1-2*tau_bar/m).*Ytr + 1);

            ub(1:n,1) = inf; %C
            
            for k = 1 : m
                ub(n*k+1:n*(k+1),1) = 0.5*((1-2*tau(k)).*Ytr + 1);
            end
            
            % equality constraints
            
            A1 = eye(n);
            
            A2Cell = repmat({-ones(1,m)}, 1, n);
            A2 = blkdiag(A2Cell{:});
            
            Aeq1 = [A1 A2];
            beq1 = zeros(n,1);
            
            A3 = zeros(m,n);
            
            A4Cell = repmat({Ytr'}, 1, m);
            A4 = blkdiag(A4Cell{:});
            
            Aeq2 = [A3 A4];
            beq2 = zeros(m,1);
            
            Aeq = [Aeq1 ; Aeq2];
            beq = [beq1 ; beq2];
            %%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXXXXXXXXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            H = sparse(H); Aeq = sparse(Aeq); beq = sparse(beq);
            
            options = optimset('Display', 'on');
            
            x = quadprog(H,-f,[],[],Aeq,beq,lb,ub,[],options,'Diagnostics','on');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXXXXXXXXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            alpha_bar = x(1:n,:);

            alpha_k = {};

            for k = 1 : m
                alpha_k{k} = x(n*k+1:n*(k+1),:);
            end
     

%%%%%%%%%%%%%%%%%%%%%%%%%%% XXXXXXXXXXXXXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% alpha1 = x(n+1:end,1);
% num_rows = numel(alpha1) / m; % Calculate the number of sets of m rows
% 
% % Reshape the column vector into a matrix with m rows
% reshaped_matrix = reshape(alpha1, m, num_rows);
% 
% % Sum along the rows (sum each set of m rows)
% sums = sum(reshaped_matrix, 1);
% 
% sum(sums==alpha_bar)

