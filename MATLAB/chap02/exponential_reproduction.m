clear all
% close all

% Font sizes
label_size = 28;
axis_size  = 24;
text_size  = 28;

% Sampling period in seconds
T = 1;

% "Continuous" time resolution
T_s = T/256;

% E-spline parameters
P       = 6; 
m       = (0:P)';
% alpha_0 = -1j * pi * P / (P + 1);
% lambda  = 1j * 2 * pi / (P + 1);
alpha_0 = -1j * pi /2;
lambda  = 1j * pi / P;
alpha_vec = alpha_0 + lambda * m;
[phi, t_phi] = generate_e_spline(alpha_vec, T_s, T); % compute in time
phi = real(phi);

% Number of samples
N = 16;

% Plot the B-spline
max_x = max(phi);
min_x = min(phi);
dyn_r = max_x - min_x;

figure
set(gcf, 'Position', [50 50 420 315])
plot(t_phi, phi, 'k', 'LineWidth', 2)
xlabel('$$t$$ [s]', 'Interpreter', 'Latex', 'FontSize', text_size)
set(gca, 'FontSize', axis_size)
axis([t_phi(1) t_phi(end) min_x-.2*dyn_r max_x+.4*dyn_r])

w = -2.8*pi:0.01:2.8*pi;
[X, Y]           = meshgrid(alpha_vec, 1j*w);
num              = 1 - exp(X - Y);
denum            = Y - X;
indet_idx        = (num == 0) & (denum == 0);
num(indet_idx)   = 1;
denum(indet_idx) = 1;
if numel(alpha_vec) > 1
    beta_w = abs(prod(num./denum, 2));
else
    beta_w = abs(prod(num./denum));
end

beta_w    = 20*log10(beta_w);
min_b     = -min(beta_w);
dyn_range = max(beta_w) - min(beta_w);
tick_0    = mod(dyn_range, 50);

figure
set(gcf, 'Position', [200 50 420 315])
plot(w, min_b+beta_w, 'k', 'LineWidth', 2)
hold on
hdl = xlabel('$$\omega$$ [rad$$\cdot$$s$$^{-1}]$$', 'Interpreter', 'Latex');
set(hdl, 'FontSize', label_size)
hdl = ylabel('[dB]');
set(hdl, 'FontSize', label_size)
ticks_y = tick_0 + (0 : 50 : 150);
set(gca, 'YTick', ticks_y);
set(gca, 'YTickLabel', {'-150', '-100', '-50', '0'})
set(gca, 'XTick', -4*pi:1*pi:4*pi);
set(gca, 'XTickLabel', {'-4\pi', '-3\pi', '-2\pi', '-\pi', '0', '\pi', '2\pi', '3\pi', '4\pi'});
set(gca, 'FontSize', axis_size)
axis([w(1) w(end) 0 max(min_b+beta_w)+.2*dyn_range])

figure
set(gcf, 'Position', [350 50 420 315])
r   = 1;
ang = 0:0.01:2*pi; 
xp  = r*cos(ang);
yp  = r*sin(ang);
plot(xp, yp, 'k')
hold on
plot(real(exp(alpha_vec)), imag(exp(alpha_vec)), '.k', 'MarkerSize', 30)
hdl = text(real(1.2*exp(alpha_vec(1))), imag(1.2*exp(alpha_vec(1))), '$$e^{\alpha_0}$$');
set(hdl, 'FontSize', text_size)
set(hdl, 'Interpreter', 'Latex');
hdl = text(real(1.2*exp(alpha_vec(2))), imag(1.2*exp(alpha_vec(2))), '$$e^{\alpha_1}$$');
set(hdl, 'FontSize', text_size)
set(hdl, 'Interpreter', 'Latex');
hdl = text(real(1.2*exp(alpha_vec(end))), imag(1.2*exp(alpha_vec(end))), '$$e^{\alpha_P}$$');
set(hdl, 'FontSize', text_size)
set(hdl, 'Interpreter', 'Latex');
hdl = xlabel('Real Part', 'Interpreter', 'Latex');
set(hdl, 'FontSize', label_size)
hdl = ylabel('Imaginary Part', 'Interpreter', 'Latex');
set(hdl, 'FontSize', label_size)
set(gca, 'FontSize', axis_size)
set(gca, 'XTick', [-1 0 1]);
set(gca, 'YTick', [-1 0 1]);
axis([-1.4 1.4 -1.4 1.4])
axis square
hdl = line([-2; 2], [0; 0]); set(hdl, 'Color', 'k', 'LineStyle', ':')
hdl = line([0; 0], [-2; 2]); set(hdl, 'Color', 'k', 'LineStyle', ':')

% Time interval that we want to reconstruct
if mod(N, 2) == 0
    n1 = -N/2;
    n2 = N/2 - 1;
else
    n1 = -(N-1)/2;
    n2 = (N-1)/2;
end
n_vec = (n1:n2)' - floor((P+1)/2);
t1    = n_vec(1) * T + t_phi(1);
t2    = n_vec(end) * T + t_phi(end);
t     = (t1:T_s:t2)';
L     = length(t);
t_int = t2 - t1 + T_s;

% c_m_n coefficients
c_m_n = get_c_m_n_exp(alpha_vec, n_vec, phi, t_phi, T);

for ith_m = 1:P+1
    m = ith_m - 1;
    
    alpha = alpha_vec(ith_m);
    
    % Original exponential
    f = exp(alpha * t / T);
    
    % Exponential reproduction
    f_rec = zeros(length(t), length(n_vec)+1);
    for i = 1:length(n_vec)
        [~, idx_f, idx_phi] = intersect(round(t/T_s), round((t_phi+n_vec(i)*T)/T_s));
        f_rec(idx_f,i+1)    = c_m_n(ith_m, i) * phi(idx_phi);

        f_rec(:,1) = f_rec(:,1) + f_rec(:,i+1);
    end

    % Dynamic range of the reproduced polynomial
    dyn_range = max(real(f_rec(:,1))) - min(real(f_rec(:,1)));
    
    % Measure the error
    phi_support = length(phi);
    idx         = (1+phi_support:length(f)-phi_support)';
    rec_er      = f(idx) - f_rec(idx,1);
    MSE         = (rec_er' * rec_er) / length(rec_er);
    disp(['MSE = ' num2str(MSE)])

    figure
    set(gcf, 'Position', [50+100*(ith_m-1) 500 420 315])
    plot(t, real(f_rec(:,1)), 'k', 'LineWidth', 4)
    hold on
    for i = 1:length(n_vec)
        plot(t, real(f_rec(:,i+1)), '-k')
    end
    plot(t, f, '--k', 'LineWidth', 2)
    hdl = xlabel('$$t$$ [s]', 'Interpreter', 'Latex');
    set(hdl, 'FontSize', label_size)
    set(gca, 'FontSize', axis_size)
    axis([t(1) t(end) min(real(f_rec(:,1)))-.3*dyn_range max(real(f_rec(:,1)))+.3*dyn_range])
    
end

