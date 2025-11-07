function PSOtutor

clear; close all; clc; 

%% Objective Functions (Pick one)
% Himmelblau Function
%     func = @(x) (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;
%     xLow = [-5 -5];             % Lower Bounds [x1_low, x2_low]
%     xHi = [5 5];                % Upper Bounds [x1_hi,  x2_hi]

% 2-D Rastrigin Function
    func = @(x) 20 + sum(x.^2) - 10*sum(cos(2*pi*x));
    xLow = [-5.12 -5.12];       % Lower Bounds [x1_low, x2_low]
    xHi = [5.12 5.12];          % Upper Bounds [x1_hi,  x2_hi]

% Peaks Function
%     func = @(x) peaks(x(1),x(2));
%     xLow = [-3 -3];             % Lower Bounds [x1_low, x2_low]
%     xHi = [3 3];                % Upper Bounds [x1_hi,  x2_hi]

% Schwefel Function
%     func = @(x) -sum(x.*sin(sqrt(abs(x))));
%     xLow = [-500 -500];           % Lower Bounds [x1_low, x2_low]
%     xHi = [500 500];              % Upper Bounds [x1_hi,  x2_hi]

%% Major Parameters [MAKE CHANGES HERE]
    tstart = tic;       % Starts the timer
    N = 2;              % No. of variables
    P = 100;            % No. of particles
    I = 100;            % Max no. of iterations (exit criterion #1)
    tol = 1e-3;         % Tolerance for convergence (exit criterion #2)
    w = [0.9 0.1];      % Linearly decreasing inertia weight, c0
    xR = xHi - xLow;    % Range of Bounds
    vM = xR*0.1;        % Maximum velocity
    c1 = 2;             % Cognitive Parameter
    c2 = 2;             % Social Parameter
 
%% Initialization
    S = repmat(struct('pos',0,...       % Position vector, x
                      'obj',0,...       % Objective function at pos
                      'bestPos',0,...   % Best position so far
                      'bestObj',0,...   % Best obj. func. so far
                      'vel',0),...      % Particle's velocity
                      P,1);
    for j = 1:P
        S(j).pos = rand(1,N).*xR + xLow;    % Set a random position
        S(j).vel = (2*rand(1,N) - 1).*vM;   % Set a random velocity
        S(j).obj = func(S(j).pos);          % Compute the obj. func.
        S(j).bestPos = S(j).pos;            % Set personal best position
        S(j).bestObj = S(j).obj;            % Set best obj. func.
    end
    
%% Swarming
    gen = cell(I,1);                              % Save all the swarms
    xG = zeros(I+1,3);                            % Save global best in xG
    xG(1,1:2) = rand(1,N).*xR + xLow;             % Set random point as xG
    xG(1,3) = func(xG(1,1:2));
    for iter = 1:I                                %  (exit criterion #1)
        f = vertcat(S.obj);
        if all(abs(f - min(f)) < tol)             %  (exit criterion #2)
            gen(iter:end) = []; 
            xG(iter:end,:) = []; break; 
        end 
        gen{iter} = S;                            % Save the current swarm
        
        c0 = w(1) - (w(1)-w(2))*(iter-1)/(I-1);   % Decrease c0
        for j = 1:P
            v = c0*S(j).vel ...                         % Inertia
                + c1*rand*(S(j).bestPos - S(j).pos)...  % Cognitive
                + c2*rand*(xG(iter,1:2) - S(j).pos);    % Social
            S(j).vel = sign(v).*min(abs(v),vM);   % Clamp to max velocity            
            S(j).pos = S(j).pos + S(j).vel;       % Move the particle
            S(j).obj = func(S(j).pos);            % Compute new obj. func
            if S(j).obj < S(j).bestObj            % Update personal best
                S(j).bestPos = S(j).pos;
                S(j).bestObj = S(j).obj;
            end
            if S(j).obj < xG(iter,3)              % Update global best
                xG(iter,3) = S(j).obj;
                xG(iter,1:2) = S(j).pos;
            end
        end
        xG(iter+1,:) = xG(iter,:);
    end
    I = length(gen);
    fprintf('Stopped at Iteration %d\n',I);
    fprintf('Minimum:\n f( %.6f , %.6f ) = %e\n',xG(end,:));
    toc(tstart);                                  % End timer
    
%% Plotting the Swarms
    close all;
    x1 = linspace(xLow(1),xHi(1));
    x2 = linspace(xLow(2),xHi(2));
    [X,Y] = meshgrid(x1,x2); Z = [X(:) Y(:)];
    Z = reshape(arrayfun(@(j) func(Z(j,:)),1:numel(X)),size(X));
    for j = 1:I
        contourf(X,Y,Z,100,'LineColor','none'); 
        colormap('jet'); hold on; colorbar; 
        pos = vertcat(gen{j}.pos);
        scatter(pos(:,1),pos(:,2),25,'m','filled',...
            'MarkerEdgeColor','k');
        scatter(xG(j,1),xG(j,2),25,'g','filled',...
            'MarkerEdgeColor','k');
        title(sprintf('Iteration: %d out of %d',j,I));
        axis([xLow(1) xHi(1) xLow(2) xHi(2)]);
        pause(0.1); if j < I, clf; end
    end
    
%% Swarm Performance
    figure(2); record = zeros(I,3); col = 'rbg';
    for i = 1:I
        record(i,1) = max([gen{i}.obj]);            % worst
        record(i,2) = median([gen{i}.obj]);         % median
        record(i,3) = min([gen{i}.obj]);            % best
    end
    neg = any(record <= 0,'all');
    for j = 1:3
        if neg, plot(1:I,record(:,j),col(j),'LineWidth',1.5); 
        else,   semilogy(1:I,record(:,j),col(j),'LineWidth',1.5);
        end
        hold on;
    end
    hold off; xlabel('Generation'); ylabel('Objective Function Value'); 
    title('Performance'); grid on;  legend('Worst','Median','Best');
end