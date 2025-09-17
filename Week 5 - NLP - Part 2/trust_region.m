% Reference: 
% https://github.com/anweshpanda/Trust_Region/blob/main/Trust%20Region.ipynb

% Rosenbrock function example
f = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
df = @(x) [200*(x(2)-x(1)^2)*-2*x(1) - 2*(1-x(1));
           200*(x(2)-x(1)^2)]; % gradient
df2 = @(x) [1200*x(1)^2-400*x(2)+2, -400*x(1);
            -400*x(1), 200];   % Hessian

[X, Y] = meshgrid(linspace(-1.5, 1.5, 200), ...
                  linspace(-0.5, 3, 200));
Z = 100*(Y-X.^2).^2 + (1-X).^2;
x0 = [Inf; Inf];
x1 = [-1.2; 1.0];

x = x1; % For collecting x per iteration
contourf(X, Y, Z, 50, 'FaceAlpha', 0.6); hold on;
scatter(x1(1), x1(2),'r','filled');

del = 1;        % Trust-region radius
max_del = 10;   % Max trust-region radius
dels = [];
ctr = 1;
while norm(x0-x1) > 1e-7
    x0 = x1;
    g = df(x0);
    H = df2(x0);
    dels(ctr) = del; %#ok
    
    % Use dogleg to solve the trust-region subproblem for p
    pN = -inv(H)*g;           % newton step
    pS = -g'*g/(g'*H*g)*g;    % steepest descent step
    if norm(pN) <= del
        p = pN;               % take the newton step
    elseif norm(pS) >= del
        p = del*pS/norm(pS);  % scale pS to trust-region radius
    else
        pd = pN - pS;
        d = (pS'*pd)^2 - pd'*pd*(pS'*pS-del^2);
        tau = (-pS'*pd + sqrt(d))/(pd'*pd) + 1;
        if tau < 1
            p = pS*tau;
        else
            p = pS + (tau-1)*pd;
        end
    end

    % Check for accuracy of approximation using rho
    rho = (f(x0) - f(x0+p))/(-g'*p + 0.5*p'*(H*p));
    if rho < 0.25             % Approx. is accurate
        del = 0.25 * del;     % Reduce trust-region radius
    else                      % Approx. is not accurate
        if rho > 0.75
            del = min(2*del, max_del); % Enlarge trust-region radius
        end
    end
    
    % Update x
    x1 = x0 + p;
    x = [x, x1];     %#ok
    ctr = ctr + 1;
end

% Plot all x and trust regions
t = linspace(0, 2*pi, 200);
for j = 1:length(dels)
    fill(dels(j)*sin(t)+x(1, j), ...
         dels(j)*cos(t)+x(2, j), 'k', ...
        'FaceAlpha', 0.2,'EdgeColor','none');
end
plot(x(1,:), x(2,:) ,'ro-');
hold off; grid on; axis equal;
axis([-1.5 1.5 -0.5 3]);
set(gcf,'Color','w');
fprintf('Minimum found: (%.4f, %.4f)\n',x1);
fprintf('No. of iterations: %d\n', ctr);