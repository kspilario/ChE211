close all;

% Illustrate the NLP
f = @(x) 1+0.1*(x(1).^2-x(2).^2);
h = @(x) x(1)^2 + x(2)^2 - 1;

[X, Y] = meshgrid(linspace(-6, 6, 1000), ...
                  linspace(-6, 6, 1000));
Z = 1+0.1*(X.^2-Y.^2); figure(1);
surf(X, Y, Z, 'FaceAlpha',0.5);
shading interp; hold on;
contour3(X, Y, Z, (-2:4)-0.1, '-y', 'LineWidth',1.2);

t = linspace(0, 2*pi, 200)';
x = 3*cos(t); y = 3*sin(t);  z = zeros(size(x));
for j = 1:length(t), z(j) = f([x(j) y(j)]); end
plot3(x, y, z, 'r', 'LineWidth',2);
colorbar; set(gcf, 'Color','w'); box on;

% Illustrate Barrier function
figure(2); t = linspace(0.01, 5, 500);
for mu = [0, 0.1, 1, 2, 5, 10]
    plot(t, -mu*log(t), 'LineWidth',1.2, ...
        'DisplayName',sprintf('\\mu = %.2f', mu)); 
    hold on;
end
grid on; legend; 
xlabel('x')
ylabel('Barrier Function $-\mu \ln(x)$', 'Interpreter','latex')
set(gcf, 'Color','w');