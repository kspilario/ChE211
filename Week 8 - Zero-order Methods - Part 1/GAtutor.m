function GAtutor

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
    P = 200;            % No. of chromosomes produced each generation
    B = 25;             % No. of bits for binary coding
    E = 10;             % No. of elites retained each generation
    G = 100;            % Max no. of generations (exit criterion #1)
    tol = 1e-3;         % Tolerance for convergence (exit criterion #2)
    Cr = 1.0;           % Crossover Rate
    Mr = 0.01;          % Mutation Rate
    x1B = [xLow(1) xHi(1)]; x2B = [xLow(2) xHi(2)];
    h1 = (x1B(2) - x1B(1))/(2^B - 1);   % interval bet. values of x1 
    h2 = (x2B(2) - x2B(1))/(2^B - 1);   % interval bet. values of x2

%% Initialization
    % Initialize chromosomes using numbers mapped to [0,2^B - 1]
    gen = zeros(G,P); 
    gen(1,:) = ceil(rand(1,P).*(2^(2*B) - 1));
    
%% (Fitness-Proportionate) Selection and Mating
    for c = 2:G                             %  (exit criterion #1)
        f = fitnessOf(gen(c-1,:));
        f = f - min(f);
        if all(abs(f) < tol)                %  (exit criterion #2)
            gen(c:end,:) = []; break;
        end
        p = f./sum(f);                      % Get selection probability
        [p,id] = sort(p,'descend');         % Sort the probabilities 
        gen(c-1,:) = gen(c-1,id);           % Apply ordering to prev gen
        p = cumsum(p);                      % Get cumulative probability
        gen(c,1:E) = gen(c-1,1:E);          % Secure E no. of elites
        m = E + 1;
        while m <= P                        % Build the next generation:
            if rand < Mr && m <= P          %  by MUTATION
                gen(c,m) = mutate(gen(c-1,find(p > rand,1)));
                m = m + 1;
            end
            if rand < Cr && m <= P          %  by CROSSOVER
                y = cross1([gen(c-1,find(p > rand,1)),...
                            gen(c-1,find(p > rand,1))]);
                gen(c,m) = y(randi(2));     % Choose a random child
                m = m + 1;
            end
        end
    end
    I = size(gen,1);
    [~,id] = sort(fitnessOf(gen(I,:)),'descend');
    gen(I,:) = gen(I,id); x = positionOf(gen(I,1));
    fprintf('Stopped at Generation %d\n',I);
    fprintf('Minimum:\n f( %.6f , %.6f ) = %e\n',x,func(x));
    toc(tstart);
    
%% Plotting the Generations
    x1 = linspace(x1B(1),x1B(2));
    x2 = linspace(x2B(1),x2B(2));
    [X,Y] = meshgrid(x1,x2); Z = [X(:) Y(:)];
    Z = reshape(arrayfun(@(j) func(Z(j,:)),1:numel(X)),size(X));
    close all; figure(1);
    for j = 1:I
        x = positionOf(gen(j,:));
        contourf(X,Y,Z,100,'LineColor','none'); 
        colormap('jet'); hold on; colorbar;
        scatter(x(2:end,1),x(2:end,2),25,'m','filled',...
            'MarkerEdgeColor','k');
        scatter(x(1,1),x(1,2),25,'g','filled',...
            'MarkerEdgeColor','k');
        title(sprintf('Iteration: %d out of %d',j,I));
        pause(0.1); if j < I, clf; end 
    end
    
%% Generation Performance
    figure(2); record = zeros(I,3); col = 'rbg';
    for j = 1:I
        record(j,1) = -fitnessOf(gen(j,P));             % worst
        record(j,2) = -fitnessOf(gen(j,floor(P/2)));    % median
        record(j,3) = -fitnessOf(gen(j,1));             % best
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
    
%% Crossover (One-Point)
    function F = cross1(y)
        t = ceil(rand*2*B);
        u = mod(y(1),2^t); v = mod(y(2),2^t);
        F = [y(1) - u + v,y(2) - v + u];
    end
%% Mutation (One-Point)
    function F = mutate(y)
        t = ceil(rand*2*B);
        F = bitset(y,t,1-bitget(y,t));
    end
%% Fitness of a Chromosome "n"
    function F = fitnessOf(n)
        F = -arrayfun(@(j) func(positionOf(j)),n);
    end
%% Position of a Chromosome "n"
    function F = positionOf(n)
        x1t = bitshift(n(:),-B); x2t = mod(n(:),2^B);
        F = [h1*x1t + x1B(1),h2*x2t + x2B(1)];
    end
end