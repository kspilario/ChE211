function [bestX, bestF, history] = simulated_annealing(f, x0)
    % Simulated Annealing for minimization
    % This is almost 100% AI-generated using the prompt:
    % "Write a MATLAB code that implements simulated annealing algorithm 
    % "from scratch. Use the Rosenbrock function minimization as test bed"
    % and
    % "In this code, the history may contain duplicate successive 
    % elements. Can you modify the code to store only the accepted history 
    % without duplicates?"

    % --- Parameters ---
    maxIter = 5000;          % number of iterations
    T0 = 1.0;                % initial temperature
    alpha = 0.995;           % cooling rate
    dim = length(x0);        % no. of dimensions
    stepSize = 0.2;          % step size for neighbor generation

    % --- Initialization ---
    x = x0; fx = f(x0);
    bestX = x;
    bestF = fx;

    % --- SA Main Loop ---
    T = T0;
    history = [x, fx];
    for k = 1:maxIter
        % Generate neighbor solution
        x_new = x + stepSize*randn(1, dim);
        f_new = f(x_new);

        % Compute acceptance probability
        dE = f_new - fx;
        if dE < 0 || rand < exp(-dE/T)
            x = x_new;
            fx = f_new;
            history(end+1, :) = [x, fx]; %#ok Store only accepted ones
        end

        % Track best solution
        if fx < bestF
            bestX = x;
            bestF = fx;
        end

        % Cooling
        T = alpha * T;
    end
end

rosen = @(x) (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2;
x0 = [-1.2, 1.0];

rng(1); % For replicability
[bestX, bestF, history] = simulated_annealing(rosen, x0);
fprintf('Best solution found: x = [%f, %f], f = %e\n', ...
        bestX(1), bestX(2), bestF);

close all;
[xg, yg] = meshgrid(linspace(-2,2,200), linspace(-1,3,200));
zg = 100*(yg - xg.^2).^2 + (1 - xg).^2;

figure;
contourf(xg, yg, zg, 30, 'FaceAlpha', 0.6); hold on;
plot(1,1,'r*','MarkerSize',12,'LineWidth',1.5); % true minimizer
colormap("parula"); xlabel('x1'); ylabel('x2');
title('Simulated Annealing on Rosenbrock');
axis tight; grid on; set(gcf, 'Color', 'w');

best_last = history(1, :);
for j = 2:size(history, 1)
    plot(history([j-1, j], 1), history([j-1, j], 2), ...
        'MarkerSize', 2, 'MarkerFaceColor','k', 'Color', [0 0 0 0.2]);
    if history(j, 3) < best_last(3)
        plot([best_last(1), history(j, 1)], ...
             [best_last(2), history(j, 2)], 'bo-', ...
             'MarkerSize', 4, 'MarkerFaceColor','b');
        best_last = history(j, :);
    end
    pause(0.1);
end
