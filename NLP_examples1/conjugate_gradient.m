% Rosenbrock function example
% f = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
% df = @(x) [200*(x(2)-x(1)^2)*-2*x(1) - 2*(1-x(1));
%            200*(x(2)-x(1)^2)];
% 
% [X, Y] = meshgrid(linspace(-1.5, 1.5, 200), ...
%                   linspace(-0.5, 3, 200));
% Z = 100*(Y-X.^2).^2 + (1-X).^2;
% x0 = [Inf; Inf];
% x1 = [-1.2; 1.0];

% Ellipse example
f = @(x) 0.06*(x(1)^2 - x(1)*x(2)) + 0.1*x(2)^2 - 0.3*(x(1)-x(2));
df = @(x) [0.06*(2*x(1)-x(2))-0.3; -0.06*x(1)+0.2*x(2)+0.3];
[X, Y] = meshgrid(linspace(-2, 3, 200), ...
                  linspace(-2, 3, 200));
Z = 0.06*(X.^2-X.*Y)+0.1*Y.^2-0.3*(X-Y);
x0 = [Inf; Inf];
x1 = [0; 0];

x = x1; % For collecting x per iteration
s = -df(x1);

contourf(X, Y, Z, 50, 'FaceAlpha', 0.6); hold on;
scatter(x1(1), x1(2),'r','filled');

ctr = 1;
while norm(x0-x1) > 1e-7

    x0 = x1; 

    % Run quasi-newton to get step size, a
    a0 = 0; a1 = 0.066; c = 1;
    dda = @(a,x,s) df(x+a*s)'*s;
    while abs(a0-a1) > 1e-7 && ~isinf(a1)
        an = a1 - (a0 - a1)*dda(a1,x0,s)/...
            (dda(a0,x0,s)-dda(a1,x0,s));
        a0 = a1; a1 = an; c = c + 1;
    end
    
    % Update x
    x1 = x0 + a0*s;

    x = [x, x1]; %#ok

    % Update search direction by conjugate gradient
    s = -df(x1) + s*(df(x1)'*df(x1))/(df(x0)'*df(x0));

    ctr = ctr + 1;
    % fprintf('(%.4f, %.4f), df = [%.4f, %.4f], a = %.4e, ctr = %d\n', ...
    %     x1, df(x1), a0, c);
end

% Plot all x
plot(x(1,:), x(2,:) ,'ro-');

hold off; grid on;
set(gcf,'Color','w');
fprintf('Minimum found: (%.4f, %.4f)\n',x1);
fprintf('No. of iterations: %d\n', ctr);