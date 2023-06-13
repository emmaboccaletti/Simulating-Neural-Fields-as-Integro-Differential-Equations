%% Cleaning
clear all, close all, clc, format short e;

%% Perform a numerical simulation of the neural field
% ∂𝑡𝑢(𝑥,𝑡) = −𝑢(𝑥,𝑡)+∫R 𝑤(𝑥−𝑦)𝑓(𝑢(𝑦,𝑡)−h)𝑑𝑦, 
% 𝑢(𝑥,0) = 𝑢0(𝑥),

% Setting up function handles 
w_handle = @(x) (1- abs(x)) .* exp(-abs(x));
f_handle = @(u) 1./(1 + exp(-10*u));
u_zero_handle = @(x) 1./(cosh(0.5 * x)).^2;

% Setting up the interval and how many x are in the interval 
L = 10;
% L = 50;
n = 3000;
Tfinal = 100;
% Tfinal = 500;
Tn = 3000;

%create n with n+1 so we have the 0
xvec = linspace(-L,L,n+1)'; 

% delete the last value of xvec (aka L) so it's a circle
xvec = xvec(1:n);
tspan = linspace(0,Tfinal,Tn+1)';
delta_x = xvec(2)-xvec(1);

% Setting up the vector U(0)
u_zero = u_zero_handle(xvec);

% Setting up the vector w(x)
w = w_handle(xvec)*delta_x;

%% Creating the matrix W
W = zeros(n,n);

% for loop
N_L0 = (n+2)/2;
W(N_L0,:) = w;
halfn = n/2;
for k = 1:(halfn)
    W((N_L0 - k), :) = circshift(w,k);
end
for k = 1:(halfn-1)
    W((N_L0 + k), :) = circshift(w,-k);
end


%% Setting up RHS with I(x,t)=0
rhs_handle1 = @(t,u) -u + W * f_handle(u);

% Check RHS with t=0
rhs_time0 = rhs_handle1(0, u_zero);

% Obtain t1 and y1
[t1,u1] = ode45(rhs_handle1,tspan,u_zero);


%% Simulate the neural field for t ∈ [0, 100]
figure(2)
[X, T] = meshgrid(xvec,t1);
surf(X, T, u1);
shading interp
title('Simulate the neural field where t∈[0, 100], S=[-10,10] and n=3000');
xlabel('x ∈ S=[-10,10]');
ylabel('Time t');
zlabel('solution u(x,t) without a stimulus');
