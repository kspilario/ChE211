close all;
[X, Y] = meshgrid(linspace(0, 7, 200), ...
                  linspace(0, 7, 200));

Z = (X-3).^2 + (Y-4).^2;
contourf(X, Y, Z, 15, 'HandleVisibility', 'off');
set(gcf, 'Color', 'w');

hold on;
fill([0 0 7 7 5],[5 7 7 0 0],'y','FaceAlpha',0.4, ...
    'HandleVisibility', 'off');
plot([0 5], [5 0], 'y', 'LineWidth', 2, 'DisplayName','5-x1-x2=0');

fill([2.5 7 7], [0 4.5 0],'m','FaceAlpha',0.4, ...
    'HandleVisibility','off');
plot([2.5 7],[0 4.5],'m','LineWidth',2, 'DisplayName','-2.5+x1-x2=0');
plot([0 0],[0 7],'r','LineWidth',2,'DisplayName','x1 = 0')
plot([0 7],[0 0],'b','LineWidth',2,'DisplayName','x2 = 0')
scatter(2, 3, 30, 'k','filled','MarkerEdgeColor','w','DisplayName','Optimal Point')
legend();