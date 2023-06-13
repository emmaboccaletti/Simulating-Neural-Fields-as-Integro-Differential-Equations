%% Cleaning
clear all, close all, clc, format short e;

%% Perform a numerical simulation of non- local integro-differential equations
% τ ∂𝑢(𝑥,𝑡)/∂t = −u(x,t) -β v(x,t) + ∫D 𝑤(𝑥−𝑦)𝑓(𝑢(𝑦,𝑡)−h)𝑑𝑦+𝐼(𝑥,𝑡), 
% 1/α 𝑢(𝑥,0)/∂t = u(x,t) -v(x,t)

% Setting up the variables for the function handles
%theta = 0.05;
theta = 0.3;
%k = 0.24;
k = 0.26;
r= 0.25;

% Setting up function handles 
lil_F = @(u) 1./(1 + exp(-r*( u - theta)));
df = @(u) r* lil_F(u) .* (1- lil_F(u));
f_handle = @(u) k * (lil_F(u) - lil_F(0))/df(0);

% Setting up the interval and how many x are in the interval 
L = 50;
n = 100;

% domain
D = linspace(-pi,pi,n+1)';
D = D(1:n);

%create n with n+1 so we have the 0
x = linspace(-L,L,n+1)'; 

% delete the last value of x (aka L) so it's a circle
x = x(1:n);

%% Information needed for ODE45 and g

% Setting up timespan
Tfinal = 50;
Tn = 3000;
tspan = linspace(0,Tfinal,Tn+1)';

% Setting up u and v for z
u_zero_handle = @(x) cos(5 * pi/L * x);
v_zero_handle = u_zero_handle;

% Setting up the vector z(0)
u_zero = u_zero_handle(x);
v_zero = v_zero_handle(x);
uv_zero = [u_zero; v_zero];

%% Setting up the w(x,y)

% Variables
%rho = 1.2;
rho = 1.05;
A = 5; 
a = 0.125; 
B = 4; 
b = 0.005;
p = 1;

% w(x,y)
w_handle = @(x) rho*(A*a^(p/2)*exp(-a*x.^2)-B*b^(p/2)*exp(-b*x.^2));

% Setting up the vector w(x)
delta_x = x(2)-x(1);
w = w_handle(x)*delta_x;

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
W;

%% Setting up RHS with I(x,t)
% This can be done in two steps, take the variables defined above and then
% a function rhs_handle which takes z,t, x, M, f, I, tau, alpha, beta
% then use it to generate a function rhs which only takes z and t

% z: vector composed by u_zero and v_zero,
% t: the timespan
% x: a location in  the ring domain D 
% M: matrix containing the the synaptic connection weight from position x' to x, w(x,y),
% f: an activation function s.t. it's output are the neurons at x 
% I: external stimuli applied to x at t
% tau: constant which fixes time units without loss of generality
% beta: strength of spike frequency adaptation
% alpha: rate of spike frequency adaptation

%% New methods using RHS
tau = 1;
%beta = 0;
beta = 0.25;
alpha = 0.1;
omega = 5;
I_variable = 6;

% c(t)
c = @(t) sin(t);

% I(x,t)
%I_handle = @(t,z) I_variable * exp(-1 * (z - c(t)).^2);
I_handle = @(t,z) 0;

% rhs which only takes t and z
rhs = @(t,z) rhs_handle(z, t, x, W, f_handle, I_handle);

% Obtain t1 and y1
[t1,y1] = ode45(rhs,tspan,uv_zero);
 
y1u = y1(:,1:n);
y1v = y1(:,n+1:2*n);

%% Simulate the nerual field for t ∈ [0, 100]: Space-time-plots-of-solutions-to-stationary-stripes
txt = ['t ∈ [0, ', num2str(Tfinal), '] with L=', num2str(L), ', n=', num2str(n)];
figure(1)
[X, T] = meshgrid(x,t1);
surf(X, T, y1u);
shading interp
title('Space-time plots of the solutions u(x,t) to standing waves ', txt);
xlabel('X');
ylabel('Time t');
zlabel('solution y for u (first n columns)');

figure(2)
[X, T] = meshgrid(x,t1);
surf(X, T, y1v);
shading interp
title('Space-time plots of the solutions v(x,t) to standing waves', txt);
xlabel('X');
ylabel('Time t');
zlabel('solution y for v (n to last columns)');


