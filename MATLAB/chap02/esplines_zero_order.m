close all
clear all

% Font sizes
label_size = 26;
axis_size  = 22;
text_size  = 26;

% Plot the FT of the E-spline
w = -7*pi:0.1:7*pi;

beta_w = abs(sin(w/2) ./ (w/2));
figure
set(gcf, 'Position', [100 50 420 315])
plot(w, beta_w, 'k', 'LineWidth', 2)
hdl = xlabel('$$\omega$$ [rad$$\cdot$$s$$^{-1}]$$', 'Interpreter', 'Latex');
set(hdl, 'FontSize', label_size)
set(gca, 'XTick', -8*pi:2*pi:8*pi);
set(gca, 'XTickLabel', {'-8\pi', '-6\pi', '-4\pi', '-2\pi', '0', ...
                        '2\pi', '4\pi', '6\pi', '8\pi'});
set(gca, 'FontSize', axis_size)
dyn_range = max(beta_w) - min(beta_w);
axis([w(1) w(end) min(beta_w) max(beta_w)+.4*dyn_range])
axp = get(gca,'Position');
annotation('line', [axp(1)+axp(3)*.5 axp(1)+axp(3)*.5],[axp(2) axp(2)+axp(4)]);

w0 = 3*pi/2;
beta_w = abs(sin((w-w0)/2) ./ ((w-w0)/2));
figure
set(gcf, 'Position', [550 50 420 315])
plot(w, beta_w, 'k', 'LineWidth', 2)
hold on
plot([0 w0], [1.1 1.1], 'x-k', 'MarkerSize', 10, 'LineWidth', 2)
hdl = text(.15*w0, 1.2, '$$\omega_0$$');
set(hdl, 'FontSize', text_size)
set(hdl, 'Interpreter', 'Latex');
hdl = xlabel('$$\omega$$ [rad$$\cdot$$s$$^{-1}]$$', 'Interpreter', 'Latex');
set(hdl, 'FontSize', label_size)
set(gca, 'XTick', -8*pi:2*pi:8*pi);
set(gca, 'XTickLabel', {'-8\pi', '-6\pi', '-4\pi', '-2\pi', '0', ...
                        '2\pi', '4\pi', '6\pi', '8\pi'});
set(gca, 'FontSize', axis_size)
dyn_range = max(beta_w) - min(beta_w);
axis([w(1) w(end) min(beta_w) max(beta_w)+.4*dyn_range])
axp = get(gca,'Position');
annotation('line', [axp(1)+axp(3)*.5 axp(1)+axp(3)*.5],[axp(2) axp(2)+axp(4)]);
