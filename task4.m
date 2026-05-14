load('phase_data.mat');

Oa = 1/(9*(2^(1/3)-1));
Ob = (2^(1/3)-1)/3;

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

text(0.38, 0.6, 'Liquid', 'FontSize', 11, 'Color', 'w');
text(3.5, 0.6, 'Vapor', 'FontSize', 11, 'Color', 'w');
text(1.2, 0.4, 'Liquid-Gas Coexistence Region', 'FontSize', 10, 'Color', 'w');
text(1.05, 1.05, 'Critical Point', 'FontSize', 10, 'Color', 'w');

xlabel('Specific Volume, V', 'FontSize', 12, 'Color', 'w');
ylabel('Pressure, P', 'FontSize', 12, 'Color', 'w');
title('P-V Phase Diagram (Redlich-Kwong EOS)', 'FontSize', 13, 'Color', 'w');

ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.TickLength = [0.02 0.05];

xlim([0.3 8]);
ylim([0 1.4]);

grid on;

T_iso = [1.2, 1.0, 0.9, 0.7];
colors = ['r', 'c', 'm', 'b'];

for i = 1:length(T_iso)
    T = T_iso(i);
    if T > 1
        V_iso = linspace(3*Ob + 0.01, 8, 500);
        P_iso = 3*T./(V_iso-3*Ob) - 9*Oa./(T^0.5.*V_iso.*(V_iso+3*Ob));
        plot(V_iso, real(P_iso), colors(i), 'LineWidth', 1.2);
        text(7.5, real(P_iso(end)), ['T=' num2str(T)], 'Color', colors(i), 'FontSize', 8);
    else
        T_idx = find(abs(T_values(:) - T) < 1e-6, 1);
        if ~isempty(T_idx)
            V1 = real(V1_results(T_idx));
            V2 = real(V2_results(T_idx));
            P_eq = real(P_results(T_idx));
            V_liq = linspace(3*Ob + 0.01, V1, 200);
            V_gas = linspace(V2, 8, 200);
            P_liq = 3*T./(V_liq-3*Ob) - 9*Oa./(T^0.5.*V_liq.*(V_liq+3*Ob));
            P_gas = 3*T./(V_gas-3*Ob) - 9*Oa./(T^0.5.*V_gas.*(V_gas+3*Ob));
            plot(V_liq, real(P_liq), colors(i), 'LineWidth', 1.2);
            plot([V1 V2], [P_eq P_eq], colors(i), 'LineWidth', 1.2);
            plot(V_gas, real(P_gas), colors(i), 'LineWidth', 1.2);
            text(7.5, real(P_gas(end)), ['T=' num2str(T)], 'Color', colors(i), 'FontSize', 8);
        end
    end
end
