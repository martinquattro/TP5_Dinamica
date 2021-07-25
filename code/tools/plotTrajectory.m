function plotTrajectory(x1,y1,x2, y2, img)
    figure;
    hold on;
    plot(x1, y1, 'Linewidth', 1.5)
    plot(x2, y2, 'Linewidth', 1.5)
    legend({'Sin termino disipativo' 'Con termino disipativo'},...
        'Interpreter','Latex', "location", "northeast");
    xlabel('Tiempo [seg]','Interpreter','Latex')
    ylabel('q(t) [rad]','Interpreter','Latex')
    grid on; grid minor;
    print(strcat('./plots/trayectoria_', num2str(img), '.png'),'-dpng');
end