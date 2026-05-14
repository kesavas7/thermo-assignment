load('phase_data.mat');

V_all = [V1_results; V2_results];
P_all = [P_results; P_results];

valid = ~isnan(V_all) & ~isnan(P_all);
V_all = V_all(valid);
P_all = P_all(valid);

[V_sorted, idx] = sort(V_all);
P_sorted = P_all(idx);

figure;
set(gcf, 'Color', 'k');
ax = gca;
set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'MinorGridColor', 'w');

plot(V_sorted, P_sorted, 'w-o', 'LineWidth', 1.5, 'MarkerSize', 4);
hold on;
plot(1, 1, 'w*', 'MarkerSize', 10, 'LineWidth', 2);

xlabel('Specific Volume, V', 'FontSize', 12, 'Color', 'w');
ylabel('Pressure, P', 'FontSize', 12, 'Color', 'w');
title('P-V Phase Diagram (Redlich-Kwong EOS)', 'FontSize', 13, 'Color', 'w');

text(0.38, 0.6, 'Liquid', 'FontSize', 11, 'Color', 'w');
text(3.5, 0.6, 'Vapor', 'FontSize', 11, 'Color', 'w');
text(1.2, 0.4, 'Liquid-Gas Coexistence Region', 'FontSize', 10, 'Color', 'w');
text(1.05, 1.05, 'Critical Point', 'FontSize', 10, 'Color', 'w');

ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.TickLength = [0.02 0.05];

xlim([0.3 8]);
ylim([0 1.2]);

grid on;
