% Rosenbrock function example
close all;
[X, Y] = meshgrid(linspace(-1.5, 1.5, 200), ...
                  linspace(-0.5, 3, 200));
Z = 100*(Y-X.^2).^2 + (1-X).^2;

subplot(221);
a = surf(X, Y, Z); shading interp; box on;
hold on; view(2); a.FaceAlpha = 0.6;

% Define candidates for grid search
[Xg, Yg] = meshgrid(linspace(-1.5, 1.5, 5), ...
                    linspace(-0.5, 3, 5));
Zg = 100*(Yg-Xg.^2).^2 + (1-Xg).^2; 
Xg = Xg(:); Yg = Yg(:); Zg = Zg(:);
scatter3(Xg, Yg, Zg, 'k', 'filled');
axis image;

% Define candidates for random search
subplot(222); rng(0);
a = surf(X, Y, Z); shading interp; box on;
hold on; view(2); a.FaceAlpha = 0.6;
Xr = rand([25, 1])*3-1.5;
Yr = rand([25, 1])*3.5-0.5;
Zr = 100*(Yr-Xr.^2).^2 + (1-Xr).^2;
scatter3(Xr(:), Yr(:), Zr(:), 'k', 'filled');
axis image;

% Plot progress for grid search and random search
subplot(2, 2, [3 4]);
semilogy(1:length(Zg), Zg, 'b', 'LineWidth',1.2, ...
    'DisplayName', 'Grid Search'); hold on;
semilogy(1:length(Zr), Zr, 'r', 'LineWidth',1.2, ...
    'DisplayName', 'Random Search'); grid on; legend;
[min_g, ind_g] = min(Zg);
[min_r, ind_r] = min(Zr);

% Report the minima found by each method
fprintf('Grid search min: %.2f at (%.2f, %.2f)\n', ...
    min_g, Xg(ind_g), Yg(ind_g));
fprintf('Random search min: %.2f at (%.2f, %.2f)\n', ...
    min_r, Xr(ind_r), Yr(ind_r));

scatter(ind_g, min_g, 'b', 'filled', 'MarkerEdgeColor','k', ...
    'DisplayName','Min. found by Grid Search');
scatter(ind_r, min_r, 'r', 'filled', 'MarkerEdgeColor','k', ...
    'DisplayName','Min. found by Random Search');
xlabel('Optimization Progress'); ylabel('Objective Function')
xticks(1:length(Zg));

subplot(221);
scatter3(Xg(ind_g), Yg(ind_g), Zg(ind_g), 'b', 'filled');
subplot(222);
scatter3(Xr(ind_r), Yr(ind_r), Zr(ind_r), 'r', 'filled');
set(gcf, 'Color', 'w');