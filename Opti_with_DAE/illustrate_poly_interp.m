
tm = [0; roots([6 -6 1]); 1]; % Roots of shifted Legendre polynomial
disp('Collocation pts.')
disp(tm); rng(0);

t = linspace(0, 1, 200); n = 10;
xi = randi(6,[n 4])-3;
xi = [xi; xi(1,:)];

% Basis functions
phi1 = @(t) (t-tm(2))./(tm(1)-tm(2)).*(t-tm(3))./(tm(1)-tm(3)).*(t-tm(4))./(tm(1)-tm(4));
phi2 = @(t) (t-tm(1))./(tm(2)-tm(1)).*(t-tm(3))./(tm(2)-tm(3)).*(t-tm(4))./(tm(2)-tm(4));
phi3 = @(t) (t-tm(1))./(tm(3)-tm(1)).*(t-tm(2))./(tm(3)-tm(2)).*(t-tm(4))./(tm(3)-tm(4));
phi4 = @(t) (t-tm(1))./(tm(4)-tm(1)).*(t-tm(2))./(tm(4)-tm(2)).*(t-tm(3))./(tm(4)-tm(3));

xi_fine = [];
for j = 1:4
    xi_fine(:,j) = interp1(0:n, xi(:,j), ...
        linspace(0,n,300), 'cubic'); %#ok
end
filename = 'poly_interp.gif'; f = figure;
for j = 1:height(xi_fine)

    % Interpolation of x(t)
    x = @(t) xi_fine(j,1)*phi1(t) + xi_fine(j,2)*phi2(t) + ...
             xi_fine(j,3)*phi3(t) + xi_fine(j,4)*phi4(t);
    
    plot(t, x(t), 'k', 'LineWidth', 1.2); hold on; set(gcf, 'Color','w'); 
    scatter(tm, xi_fine(j,:), 'm', 'filled'); hold off; grid on;
    title(sprintf('x = [%.2f, %.2f, %.2f, %.2f]',xi_fine(j,:)))
    ylim([-5 5]); pause(0.2);
    exportgraphics(f, filename, 'Append', true);
end
