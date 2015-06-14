clear all
close all

% Font sizes
label_size = 28;
axis_size  = 24;
text_size  = 28;

T_s = 1/2^6;
N   = 2^14;
w   = 2*pi/(N*T_s) * (-(N/2) : (N/2 - 1))';
t   = (0 : (N-1))' * T_s;

b0_w = (1 - exp(-1i * w)) ./ (1i * w);
b0_w(w==0) = 1;
b0_t = zeros(size(t));
b0_t(t>=0 & t < 1) = 1;

b1_w = b0_w.^2;
b1_t = (1/T_s) * ifft([b1_w(end/2+1:end); b1_w(1:end/2)]);

b2_w = b0_w.^3;
b2_t = (1/T_s) * ifft([b2_w(end/2+1:end); b2_w(1:end/2)]);

max_x = max([max(real(b0_w)) max(imag(b0_w))]);
min_x = min([min(real(b0_w)) min(imag(b0_w))]);
dyn_r = max_x - min_x;

figure
set(gcf, 'Position', [50 50 420 315])
plot(w, real(b0_w), 'b', w, imag(b0_w), 'r', 'LineWidth', 2)
xlabel('$$\omega$$ [rad$$\cdot$$s$$^{-1}$$]', 'Interpreter', 'Latex', 'FontSize', text_size)
set(gca, 'XTick', -8*pi:2*pi:8*pi);
set(gca, 'XTickLabel', {'-8\pi', '-6\pi', '-4\pi', '-2\pi', '0', '2\pi', '4\pi', '6\pi', '8\pi'});
set(gca, 'FontSize', axis_size)
axis([-7*pi 7*pi min_x-.2*dyn_r max_x+.4*dyn_r])

figure
set(gcf, 'Position', [550 50 420 315])
plot([-1; -T_s; t], real([0; 0; b0_t]), 'k', 'LineWidth', 2)
xlabel('$$t$$ [s]', 'Interpreter', 'Latex', 'FontSize', text_size)
set(gca, 'FontSize', axis_size)
axis([-1 5 -.1 1.1])

max_x = max([max(real(b1_w)) max(imag(b1_w))]);
min_x = min([min(real(b1_w)) min(imag(b1_w))]);
dyn_r = max_x - min_x;

figure
set(gcf, 'Position', [50 150 420 315])
plot(w, real(b1_w), 'b', w, imag(b1_w), 'r', 'LineWidth', 2)
xlabel('$$\omega$$ [rad$$\cdot$$s$$^{-1}$$]', 'Interpreter', 'Latex', 'FontSize', text_size)
set(gca, 'XTick', -8*pi:2*pi:8*pi);
set(gca, 'XTickLabel', {'-8\pi', '-6\pi', '-4\pi', '-2\pi', '0', '2\pi', '4\pi', '6\pi', '8\pi'});
set(gca, 'FontSize', axis_size)
axis([-7*pi 7*pi min_x-.2*dyn_r max_x+.4*dyn_r])

figure
set(gcf, 'Position', [550 150 420 315])
plot([-1; -T_s; t], real([0; 0; b1_t]), 'k', 'LineWidth', 2)
xlabel('$$t$$ [s]', 'Interpreter', 'Latex', 'FontSize', text_size)
set(gca, 'FontSize', axis_size)
axis([-1 5 -.1 1.1])

max_x = max([max(real(b2_w)) max(imag(b2_w))]);
min_x = min([min(real(b2_w)) min(imag(b2_w))]);
dyn_r = max_x - min_x;

figure
set(gcf, 'Position', [50 250 420 315])
plot(w, real(b2_w), 'b', w, imag(b2_w), 'r', 'LineWidth', 2)
xlabel('$$\omega$$ [rad$$\cdot$$s$$^{-1}$$]', 'Interpreter', 'Latex', 'FontSize', text_size)
set(gca, 'XTick', -8*pi:2*pi:8*pi);
set(gca, 'XTickLabel', {'-8\pi', '-6\pi', '-4\pi', '-2\pi', '0', '2\pi', '4\pi', '6\pi', '8\pi'});
set(gca, 'FontSize', axis_size)
axis([-7*pi 7*pi min_x-.2*dyn_r max_x+.4*dyn_r])

figure
set(gcf, 'Position', [550 250 420 315])
plot([-1; -T_s; t], real([0; 0; b2_t]), 'k', 'LineWidth', 2)
xlabel('$$t$$ [s]', 'Interpreter', 'Latex', 'FontSize', text_size)
set(gca, 'FontSize', axis_size)
axis([-1 5 -.1 1.1])
