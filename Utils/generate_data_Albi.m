function [X,Y] = generate_data_Albi(N_points, p_A, p_O, mu1, mu2, Sigma1, Sigma2)

p_B = 1-p_A;

X_red = []; Y_red = [];
X_blue = []; Y_blue = [];

for i = 1:N_points

   if rand(1) < p_A

       if rand(1) < p_O % outlier

            [x_i, y_i] = MixGaussAlbi(mu1,Sigma1,1);
            y_i(y_i==0) = -1;

       else

            [x_i, y_i] = MixGaussAlbi(mu1,Sigma1,1);
            y_i(y_i==0) = +1;

       end

       X_red = [X_red; x_i];
       Y_red = [Y_red; y_i];

   else

       if rand(1) < p_O % outlier

            [x_i, y_i] = MixGaussAlbi(mu2,Sigma2,1);
            y_i(y_i==0) = +1;

       else

            [x_i, y_i] = MixGaussAlbi(mu2,Sigma2,1);
            y_i(y_i==0) = -1;

       end

       X_blue = [X_blue; x_i];
       Y_blue = [Y_blue; y_i];

   end

end

X = [X_red; X_blue];
Y = [Y_red; Y_blue];


