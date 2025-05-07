syms y(x)
eqn = diff(y) == y*x;
cond = y(1) == 1;
ySol(x) = dsolve(eqn,cond)

x = linspace(0,1,20);
plot(x, ySol(x))


