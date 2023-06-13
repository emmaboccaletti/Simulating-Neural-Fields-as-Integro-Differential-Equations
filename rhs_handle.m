function [g] = rhs_handle(z,t, x, M, f, I, tau, alpha, beta)

arguments
    z double
    t double
    x {mustBeVector}
    M double
    f {mustBeA(f, 'function_handle')}
    I {mustBeA(I, 'function_handle')}
    tau {mustBeNonempty} = 1
    alpha {mustBeNonempty} = 1
    beta {mustBeNonempty} = 1
end

[numRows,numCols] = size(z);
N = numRows/2;
u = z(1:N);
v = z(N+1:numRows);


top_g = (-u - beta * v + M * f(u) + I(t,x))/tau ;
bottom_g = (u - v)*alpha;

g = [top_g; bottom_g];

end