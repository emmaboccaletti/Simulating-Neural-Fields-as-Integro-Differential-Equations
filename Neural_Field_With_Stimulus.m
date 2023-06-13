%% Cleaning
clear all, close all, clc, format short e;

%% Perform a numerical simulation of the neural field
% ∂𝑡𝑢(𝑥,𝑡) = −𝑢(𝑥,𝑡)+∫R 𝑤(𝑥−𝑦)𝑓(𝑢(𝑦,𝑡)−h)𝑑𝑦+𝐼(𝑥,𝑡), 
% 𝑢(𝑥,0) = 𝑢0(𝑥),

% Setting up function handles 
w_handle = @(x) (1- abs(x)) .* exp(-abs(x));
f_handle = @(u) 1./(1 + exp(-10*u));
u_zero_handle = @(x) 1./(cosh(0.5 * x)).^2;

% Setting up the interval and how many x are in the interval 
L = 10;
n = 3000;
Tfinal = 100;
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

%% Setting up RHS with I(x,t)
omega = 5;
I_variable = 3;

% c(t)
c = @(t) 2*(sin((2*pi / omega)*t) + 1/3*sin((6*pi / omega)*t) + 1/5*sin(((10*pi / omega)*t))) ;

% I(x,t)
I_handle = @(t,z) I_variable * exp(-1 * (z - c(t)).^2);

% RHS
rhs_handle_withI = @(t,u) -u + W * f_handle(u) + I_handle(t,xvec);


% Check c(t) and I(t,x) with t=0
c0 = c(0);
I0 = I_handle(0,xvec);

figure(1)
plot(c(tspan), tspan);
title('Plot of c as time progresses t ∈ [0, 100]');
ylabel('Time t');
xlabel('Solution c(t)');

figure(2)
[X1, T1] = meshgrid(xvec,tspan);
z = I_handle(T1,X1);
surf(X1, T1, z);
shading interp
title('Plotting the stimulus over time and space');
xlabel('x ∈ S=[-L,L]');
ylabel('Time t ∈ [0, 100]');
zlabel('Solution of I(x,t)');


% Check RHS with t=0
rhs_withI_time0 = rhs_handle_withI(0, u_zero);

% Obtain t2 and y2
[t2,y2] = ode45(rhs_handle_withI,tspan,u_zero);

%% Simulate the neural field for 𝑡 ∈ [0, 100]
figure(3)
[X, T] = meshgrid(xvec,t2);
surf(X, T, y2);
shading interp
title('Simulate the neural field where t∈[0, 100], S=[-L,L] and n=3000');
xlabel('x ∈ S=[-L,L]');
ylabel('Time t');
zlabel('Solution u(x,t) with stimulus');

