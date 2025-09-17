[X, Y] = meshgrid(linspace(-3, 3, 200), ...
                  linspace(-2, 7, 300));
Z = @(X, Y) 4 + 4.5*X - 4*Y + X.^2 + 2*Y.^2 - 2*X.*Y + X.^4 - 2*Y.*X.^2;
contour(X, Y, Z(X, Y), 50);
hold on;
scatter([1.941, -1.053, 0.6117],[3.854, 1.028, 1.4929], ...
    20,'r','filled','MarkerEdgeColor','k')
H = @(x) [2+12*x(1)^2-4*x(2), -2-4*x(1); -2-4*x(1), 4];
disp('Hessians')
disp(H([1.941, 3.854]))
disp(H([-1.053, 1.028]))
disp(H([0.6117, 1.4929]))

disp('Hessian Eigenvalues')
disp(eig(H([1.941, 3.854])))
disp(eig(H([-1.053, 1.028])))
disp(eig(H([0.6117, 1.4929])))

disp('f(x)')
disp(Z(1.941, 3.854))
disp(Z(-1.053, 1.028))
disp(Z(0.6117, 1.4929))
set(gcf, 'Color','w');