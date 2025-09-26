% Booth function example
f = @(x) (x(1)+2*x(2)-7)^2 + (2*x(1)+x(2)-5)^2;
df = @(x) [2*(x(1)+2*x(2)-7)+2*(2*x(1)+x(2)-5)*2;
           2*(x(1)+2*x(2)-7)*2+2*(2*x(1)+x(2)-5)];

[X, Y] = meshgrid(linspace(-5, 5, 500), ...
                  linspace(-5, 5, 500));
Z = zeros(size(X));
for j = 1:length(X)
    for k = 1:length(Y)
        Z(j,k) = f([X(j, k), Y(j, k)]);
    end
end
contourf(X, Y, Z, 100);
x0 = [Inf; Inf];
x1 = [-4; -4];


x = x1; % For collecting x per iteration

hold on; scatter(x1(1), x1(2),'r','filled');

ctr = 1; J = eye(length(x1)); I = J;
while norm(x0-x1) > 1e-7

    % Calculate gradient
    x0 = x1; s = -J*df(x0);

    % Run quasi-newton to get step size, a
    a0 = 0; a1 = 0.001; c = 1;
    dda = @(a,x,s) df(x+a*s)'*s;
    while abs(a0-a1) > 1e-7 && ~isinf(a1)
        an = a1 - (a0 - a1)*dda(a1,x0,s)/...
            (dda(a0,x0,s)-dda(a1,x0,s));    
        a0 = a1; a1 = an; c = c + 1;
    end

    % Update x
    x1 = x0 + a0*s; 

    % Update inverse Hessian, J
    y = df(x1) - df(x0);
    J = (I-s*y'/(y'*s))*J*(I-y*s'/(y'*s)) + s*s'/(y'*s);

    x = [x, x1]; %#ok
    ctr = ctr + 1;

    % fprintf('(%.4f, %.4f), df = [%.4f, %.4f], a = %.4f, ctr = %d\n', ...
    %     x1, df(x1), a0, c);
end

% Plot all x
plot(x(1,:), x(2,:) ,'ro-');

hold off; grid on;
set(gcf,'Color','w');
fprintf('Minimum found: (%.4f, %.4f)\n',x1);
fprintf('No. of iterations: %d\n', ctr);