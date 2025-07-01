function [dydt] = ode_wolbachia(t, y,  b, d, d_prime, phi, h, k, C, u, v)
% ODE definition of wildtype/wolbachia system

%disp([b, d, d_prime, phi, h, k, C, u, v])

N_m = y(1);
N_w = y(2);

% X is the total population size
X = N_m + N_w;

% Y_m is the maximum number of wildtype offspring
Y_m = 0;
if X > 0
    Y_m = (N_m * (N_m + (1-v) * phi * N_w) + N_w * ((1-u) * N_m + (1-v) * phi * N_w)) / (N_m + N_w);
end

% Z_w is the maximum number of wolbachia-infected offspring
Z_w = v * phi * N_w;

dydt = [
    b*Y_m*F(X,h,k,C) - d*N_m;
    b*Z_w*F(X,h,k,C) - d_prime*N_w;
    ];
end
    
% function to model density dependent reproduction (Dye exponential form)
% X is total population size, C is maximum population size
function out = F(X, h, k, C)
    if X < C
        if  X > 0
            out = exp(-h*X^k);  % exponential decrease in reproduction rate
        else
           out =  1.0;
        end
    else
        out = 0; % if maximum population size is reached, reproduction rate is zero
    end
end

