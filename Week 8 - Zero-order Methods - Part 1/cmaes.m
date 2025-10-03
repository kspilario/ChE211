% This code is almost 100% AI-generated using the following prompts:
%   "Write a MATLAB code that implements CMAES from scratch using 
%    Rosenbrock function as test bed."
%   "Modify the code to stick with n=2, then plot the successive 2D 
%    gaussians every iteration."
%   "This code makes the gaussians explode towards the end of 
%    the iterations. How to avoid this?"
% 
% CMA-ES for Rosenbrock with fading Gaussian ellipses (stabilized)

clear; close all; 
rng(42); % For replicability

%% Problem setup
n = 2;                  % fixed dimension
maxGen = 100;           % max generations
sigma = 0.3;            % initial step-size
m = [-1.2; 1.0];        % initial mean

% CMA-ES parameters
lambda = 10;            % population size
mu = floor(lambda/2);   % parents
weights = log(mu+0.5) - log(1:mu)'; 
weights = weights / sum(weights);
mu_eff = 1 / sum(weights.^2);

% Adaptation params
c_sigma = (mu_eff+2)/(n+mu_eff+5);
d_sigma = 1 + 2*max(0, sqrt((mu_eff-1)/(n+1))) + c_sigma;
c_c = 4/(n+4);
c1 = 2/((n+1.3)^2 + mu_eff);
cmu = min(1-c1, 2*(mu_eff-2+1/mu_eff)/((n+2)^2 + mu_eff));
chi_n = sqrt(n)*(1-1/(4*n)+1/(21*n^2));

% Initialization
pc = zeros(n,1);
ps = zeros(n,1);
B = eye(n);
D = ones(n,1);
C = B*diag(D.^2)*B';
eigEvalIter = 0;
eigEvalFrequency = 1/(c1+cmu)/n/10;

% Rosenbrock
rosen = @(x) (1-x(1))^2 + 100*(x(2)-x(1)^2)^2;

%% Visualization setup
figure;
[xg, yg] = meshgrid(linspace(-2,2,200), linspace(-1,3,200));
zg = 100*(yg - xg.^2).^2 + (1 - xg).^2;
contourf(xg, yg, zg, 30, 'FaceAlpha', 0.6); hold on;
plot(1,1,'r*','MarkerSize',12,'LineWidth',1.5); % true minimizer
colormap("parula"); xlabel('x1'); ylabel('x2');
title('CMA-ES on Rosenbrock');
axis tight; grid on; set(gcf, 'Color', 'w');

colormap(jet);
cmap = colormap;
numColors = size(cmap,1);

%% Main loop
for gen = 1:maxGen
    % Sample offspring
    arz = randn(n, lambda);
    ary = B * (D .* arz);
    arx = m + sigma * ary;
    
    % Evaluate
    fitness = arrayfun(@(k) rosen(arx(:,k)), 1:lambda);
    [~, idx] = sort(fitness);
    
    xsel = arx(:,idx(1:mu));
    ysel = ary(:,idx(1:mu));
    zsel = arz(:,idx(1:mu));
    
    % Recombination
    y_w = ysel*weights;
    m = m + sigma*y_w;
    
    % Update evolution paths
    z_w = zsel*weights;
    ps = (1-c_sigma)*ps + sqrt(c_sigma*(2-c_sigma)...
        *mu_eff)*(B*((1./D).*z_w));
    h_sigma = norm(ps)/sqrt(1-(1-c_sigma)^(2*gen))/chi_n < (1.4+2/(n+1));
    pc = (1-c_c)*pc + h_sigma*sqrt(c_c*(2-c_c)*mu_eff)*y_w;
    
    % Covariance update
    C = (1-c1-cmu)*C + c1*(pc*pc' + (1-h_sigma)*c_c*(2-c_c)*C) ...
        + cmu * ysel * diag(weights) * ysel';
    
    % Step-size adaptation
    sigma = sigma * exp((c_sigma/d_sigma)*(norm(ps)/chi_n - 1));
    
    % --- Safeguards ---
    sigma = max(min(sigma, 1), 1e-10);   % clamp sigma
    
    % Eigen decomposition occasionally
    if gen - eigEvalIter > max(1, floor(eigEvalFrequency))
        eigEvalIter = gen;
        C = triu(C) + triu(C,1)'; % symmetrize
        [Btmp, Dtmp] = eig(C);
        Dvals = diag(Dtmp);
        % clamp eigenvalues to avoid exploding/vanishing axes
        Dvals = max(min(Dvals, 1e2), 1e-10);  
        B = Btmp;
        D = sqrt(Dvals);
    end
    
    %% Draw fading ellipse
    [V,E] = eig(C);
    t = linspace(0,2*pi,100);
    circle = [cos(t); sin(t)];
    %ell = m + sigma * V*sqrt(E)*circle;
    
    % map iteration to color index
    colorIdx = round((gen/maxGen)*(numColors-1))+1;
    col = cmap(colorIdx,:); 
    
    % fade factor (start pale, end vivid)
    fadeFactor = 0.2 + 0.8*(gen/maxGen);
    fadedCol = (1-fadeFactor)*[1 1 1] + fadeFactor*col;
    
    % cap maximum ellipse radius for visualization
    maxAxisLen = 2;
    ell = m + min(1,maxAxisLen/norm(V*sqrt(E)))*sigma*V*sqrt(E)*circle;
    
    fill(ell(1,:), ell(2,:), fadedCol, 'FaceAlpha', 0.2);
    plot(m(1), m(2), 'ko','MarkerFaceColor','k','MarkerSize',4);
    
    if mod(gen,10)==0
        drawnow;
        fprintf('Gen %d: best=%.3e, sigma=%.3e\n', ...
            gen, min(fitness), sigma);
    end
end
