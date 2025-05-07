% y` = x + sin(y/sqr(7))
% y(0.5) = 0.6
% [0.5 ,1.5]


f = @(x, y) x+sin(y/sqrt(7));
tspan = [0.5, 1.5];
y0 = 0.6;

[T, Y] = ode23(f, tspan, y0)
[t, y] = ode45(f, tspan, y0)

plot(T, Y, t, y)
grid on;