% Data
expected_direct = [0 0 0 0 0 0 25 0 0 -32 -25]
resolved_direct = [-5 -4 3 -6 15 4 28 -5 -5 -23 -23]

expected_ref = [-70 -34 54 63 45 45 54 38 -43 51 -38]
resolved_ref = [-68 -37 60 60 48 60 60 39 -36 60 -30]


% Computation
err_direct = abs(resolved_direct - expected_direct)
err_ref = abs(resolved_ref - expected_ref)
err = [err_direct err_ref]

set(0, 'DefaultLineLineWidth', 2);
cdfplot(err);
% set(gca,'LineWidth', 1.8);

hold on
plot(xlim, [0.5 0.5], '--')
% cdfplot(err_direct);
% cdfplot(err_ref);
% legend("Direct path", "Reflected path", "All included", 'Location','southeast')
xlabel("Error (in °)")
ylabel("Frequency")
title("CDF of AoA resolution accuracy")
hold off
set(0, 'DefaultLineLineWidth', 1);