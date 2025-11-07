f = @(x) x.^4-x+1;
df = @(x) 4*x.^3-1;
ddf = @(x) 12*x.^2;

% Use Newton's method
subplot(121);
x = linspace(-3, 3, 200);
plot(x, f(x), 'k', 'LineWidth',1.5); grid on;
set(gcf, 'Color', 'w'); hold on;
x0 = 3; x1 = Inf; ctr = 1;
scatter(x0, f(x0), 20, 'b','filled');
axis([-3 3 0 90]);
while abs(x0-x1) > 1e-7
    x1 = x0;
    x0 = x1 - df(x1)/ddf(x1);
    ctr = ctr + 1;
    scatter(x1, f(x1), 20, 'b','filled');
end
disp(x0);
fprintf('No. of iterations: %d\n', ctr);
title('Newton''s Method')

% Use quasi-newton method
subplot(122);
x = linspace(-3, 3, 200);
plot(x, f(x), 'k', 'LineWidth',1.5); 
grid on; hold on; axis([-3 3 0 90]);
x0 = 3; x1 = 2; ctr = 1;
scatter(x0, f(x0), 20, 'b','filled');
while abs(x0-x1) > 1e-7
    xn = x1 - (x0 - x1)*df(x1)/(df(x0)-df(x1));    
    x0 = x1; x1 = xn;
    scatter(x1, f(x1), 20, 'b','filled');
    ctr = ctr + 1;
end
disp(x0);
title('Quasi-newton Method')
fprintf('No. of iterations: %d\n', ctr);
