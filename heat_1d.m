% Date June 18, 2025

clc;
clear;
close all;

dt = 1.0e-2;
N = 100;
t_final = 0.05;
dx = 1.0 / N ;
x_arr = linspace(0, 1, N+1);
N_step = ceil(t_final/dt);
mu = dt / (dx * dx);

mat_A = zeros(N+1, N+1);
vec_u = zeros(N+1, 1);
vec_b = zeros(N+1, 1);

% Matrix assemble
for ii = 1:N+1
    if ii == 1
        mat_A(ii,ii) = 1.0;
    elseif ii < N+1
        mat_A(ii, ii-1) = -mu;
        mat_A(ii, ii)   = 1 + 2.0 * mu;
        mat_A(ii, ii+1) = -mu;
    else
        mat_A(ii, ii-1) = -2.0 * mu;
        mat_A(ii, ii) = 1 + 2.0 * mu;
        % mat_A(ii, ii) = 1.0;
    end
end

% Initialize temperature vector
for ii = 1:N+1
    xx = (ii-1) * dx;
    vec_u(ii) = g_value(xx);
end

% Time advancing
for i_time = 1:N_step
    tt = (i_time) * dt;
    for ii = 1:N+1
        if ii == 1
            vec_b(ii) = 0.0;
        elseif ii < N+1
            vec_b(ii) = vec_u(ii);
        else
            hh = 2.0 * pi * exp(-4.0*pi*pi*tt);
            vec_b(ii) = vec_u(ii) + 2.0 * mu * dx * hh;
        end
    end
    vec_u = linsolve(mat_A, vec_b);
end

plot(x_arr, vec_u);

%% function
function gval = g_value(xx)
    gval = sin(2.0*pi*xx);
end