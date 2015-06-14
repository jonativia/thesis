clear all
close all

label_size = 28;
axis_size  = 24;
text_size  = 28;

% Parameters
T     = 1;
t_int = T * 8;                % Temporal interval where Diracs can be located
K     = 4;                    % Number of Diracs
P     = 2*K - 1;              % E-spline order
N     = floor(t_int / T) + P; % number of temporal samples
T_s   = T / 64;               % "continuous" time resolution

% Stream of Diracs params
t_sig = (0 : T_s : t_int)';
itk = [150 190 350 420];
t_k = t_sig(itk).';
a_k = [1.5 1 2 3];

% Generate the continuous-time signal x(t)
x = zeros(size(t_sig));
x(itk(1:K)) = a_k(1:K);

% E-spline
m = 0:P;
omega_0 = -pi /2;
lambda  = pi / P;
% omega_0 = -pi * P / (P+1);
% lambda  = 2 * pi / (P+1);
omega_m   = omega_0 + lambda * m;
alpha_vec = 1j * omega_m;
[phi, t_phi] = generate_e_spline(alpha_vec, T_s, T, 'anticausal');
L_phi        = length(phi);
h            = real(phi(end:-1:1));
t_h          = -t_phi(end:-1:1);

% Plot sampling kernel and exp(alpha_vec)
figure
set(gcf, 'Position', [50 350 420 315])
plot(t_h, h, 'k', 'LineWidth', 2)
hdl = xlabel('$$t$$ [s]', 'Interpreter', 'Latex');
set(hdl, 'FontSize', label_size)
set(gca, 'FontSize', axis_size)
dyn_range = max(h) - min(h);
axis([t_h(1) t_h(end) min(h) max(h)+.4*dyn_range])

figure
set(gcf, 'Position', [500 350 420 315])
r   = 1;
ang = 0:0.01:2*pi; 
xp  = r*cos(ang);
yp  = r*sin(ang);
plot(xp, yp, 'k')
hold on
plot(real(exp(alpha_vec)), imag(exp(alpha_vec)), '.k', 'MarkerSize', 40)
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
hdl = line([-2; 2], [0; 0]); set(hdl, 'Color', 'k')
hdl = line([0; 0], [-2; 2]); set(hdl, 'Color', 'k')
set(gca, 'XTick', [-1 0 1]);
set(gca, 'YTick', [-1 0 1]);
axis([-1.4 1.4 -1.4 1.4])
axis square

% c_m_n parameters
n_vec = 1 : N;
t_n   = n_vec * T;
c_m_n = get_c_m_n_exp(alpha_vec, n_vec, phi, t_phi, T);

% Compute the continuous time signal y(t) = x(t) * h(t) and samples y[n]
y = conv(x, h);
t_0 = t_sig(1) + t_h(1);
t_f = t_sig(end) + t_h(end);
t_y = t_0 : T_s : t_f;
[~,idx] = intersect(t_y, t_n);
y_n = y(idx);

% Compute the moments of the signal
s_m = c_m_n * y_n;

% Locate the diracs
S = toeplitz(s_m(K+1:end), s_m(K+1:-1:1));
[~,~,V] = svd(S);
h       = V(:,end);
uu_k    = roots(h);
pp_k    = mod(angle(uu_k), 2*pi);
tt_k    = sort(T * pp_k / lambda);
phi_mat = get_phi_tk_n_mat(phi, t_phi, tt_k, n_vec, T, T_s);
aa_k    = real(phi_mat \ y_n);

disp('%%--Stream of deltas--%%');
disp(' ');
disp(['-> Original tk   : ', mat2str(t_k(1:K))]);
disp(['-> Estimated tk  : ', mat2str(tt_k)]);
disp(['-> Squared Error : ', mat2str((t_k(:) - tt_k(:)).^2)]);
disp(' ');

disp(['-> Original ak   : ', mat2str(a_k(1:K))]);
disp(['-> Estimated ak  : ', mat2str(aa_k)]);
disp(['-> Squared Error : ', mat2str((a_k(:) - aa_k(:)).^2)]);

% Fourier transform of signal x
omega = linspace(-2.4*pi/T, 2.4*pi/T, 1024)';
X_ft  = abs(exp(-1j * omega * t_k(:)') * a_k(:));

figure
set(gcf, 'Position', [50 200 420 315])
stem(t_k(1:K), a_k(1:K), '^k', 'fill', 'LineWidth', 1, 'MarkerSize', 10)
axis([t_y(1) t_y(end) -.5 4])
hdl = xlabel('$$t$$ [s]', 'Interpreter', 'Latex');
set(hdl, 'FontSize', label_size)
set(gca, 'FontSize', axis_size)

figure
set(gcf, 'Position', [500 200 420 315])
plot(t_y, y, 'k', 'LineWidth', 2)
hold on
stem(t_n, y_n, '.k', 'MarkerSize', 40)
hdl = legend(' $$y(t)$$ ', ' $$y[n]$$ ', 'Location', 'NorthEast');
set(hdl, 'Interpreter', 'Latex', 'FontSize', label_size)
hdl = xlabel('$$t$$ [s]', 'Interpreter', 'Latex');
set(hdl, 'FontSize', label_size)
set(gca, 'FontSize', axis_size)
dyn_range = max(y) - min(y);
axis([t_y(1) t_y(end) min(y) max(y)+.5*dyn_range])

figure
set(gcf, 'Position', [50 50 420 315])
plot(omega, X_ft, 'k', 'LineWidth', 2)
hold on
stem(-omega_m/T, abs(s_m), '.k', 'MarkerSize', 40)
hdl = legend('$$|\hat{x}(\omega)|$$ ', '   $$|s[m]|$$', 'Location', 'NorthEast');
set(hdl, 'Interpreter', 'Latex', 'FontSize', label_size)
hdl = xlabel('$$\omega$$ [rad$$\cdot$$s$$^{-1}]$$', 'Interpreter', 'Latex');
set(hdl, 'FontSize', label_size)
set(gca, 'XTick', [-2*pi/T -pi/T 0 pi/T 2*pi/T]);
set(gca, 'XTickLabel', {'-2\pi/T', '-\pi/T', '0', '\pi/T', '2\pi/T'});
set(gca, 'FontSize', axis_size)
dyn_range = max(X_ft) - min(X_ft);
axis([omega(1) omega(end) min(X_ft) max(X_ft)+.5*dyn_range])

figure
set(gcf, 'Position', [500 50 420 315])
stem(tt_k(1:K), aa_k(1:K), '^k', 'fill', 'LineWidth', 1, 'MarkerSize', 10)
axis([t_y(1) t_y(end) -.5 4])
hdl = xlabel('$$t$$ [s]', 'Interpreter', 'Latex');
set(hdl, 'FontSize', label_size)
% hdl = title('Reconstructed stream of Diracs', 'Interpreter', 'Latex');
% set(hdl, 'FontSize', label_size)
set(gca, 'FontSize', axis_size)
