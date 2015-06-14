clear all
close all
clc

label_size = 28;
axis_size  = 24;
text_size  = 28;

% Parameters
P          = 10;        % E-spline order
N          = 2*(P+1);   % Temporal samples
T          = 1 / N;     % sampling period
T_s        = T / 64;    % "continuous" time resolution
SNR_min    = 0;         % SNR range in dB
SNR_max    = 30;
it_per_snr = 100;       % Number of repetitions per SNR
cad_it     = 5;

% Time interval that we want to consider
n1 = 1;
n2 = N;
n_vec = (n1:n2)';
t1    = n1 * T;
t2    = (n2+1) * T - T_s;
t     = (t1:T_s:t2)';
L_t   = length(t);
t_int = t2 - t1 + T;

% Generate the stream of Diracs
itk = [150 600];
t_k = t(itk).';
a_k = [1 1.5];
K   = length(t_k);

% Generate the x(t) signal
x = zeros(size(t));
x(itk(1:K)) = a_k(1:K);

% E-spline
m = 0:P;
omega_0 = -pi * P / N;
lambda  = 2 * pi / N;
% omega_0 = -pi * P / (P+1);
% lambda  = 2 * pi / (P+1);
omega_m      = omega_0 + lambda * m;
alpha_vec    = 1j * omega_m;
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

% Compute the perfect exponential reproduction interval
t_start = t_phi(end) + n_vec(1) * T - T;
t_end   = t_phi(1)   + n_vec(end) * T + T;
disp('%%--Exponential perfect reconstruction time interval--%%');
disp(['From ' num2str(t_start) ' to ' num2str(t_end)])
disp(' ')

% c_m_n parameters
c_m_n = get_c_m_n_exp(alpha_vec, n_vec, phi, t_phi, T);

% Compute y_n as y_n = <x(t),phi(t/T-n)>
y_n = zeros(length(n_vec), 1);
for it = 1 : N
    idx_1 = round((t_phi(1) + n_vec(it) * T - t(1)) / T_s) + 1;
    idx_2 = round((t_phi(end) + n_vec(it) * T - t(1)) / T_s) + 1;
    idx   = (idx_1:idx_2)';
    
    [~,idx_x,idx_phi] = intersect(1:L_t, idx);
    y_n(it) = x(idx_x).' * phi(idx_phi);
end

P_y = y_n' * y_n / N;

% Variables to control the number of iterations
SNRs     = SNR_min : SNR_max;
SNR_vec  = repmat(SNRs, it_per_snr, 1);
SNR_vec  = SNR_vec(:)';
num_iter = length(SNR_vec);
ttk_pen  = zeros(num_iter, K);
ttk_cad  = zeros(num_iter, K);
SNRm_vec = zeros(num_iter, 1);

% Standard deviation of the noise for each realisation and noise matrix
e_mat  = randn(length(y_n), num_iter);
P_e    = sum(e_mat .* conj(e_mat)) / N;
sigmas = sqrt(10.^(-SNR_vec/10) .* repmat(P_y, 1, num_iter) ./ P_e);
e_mat  = repmat(sigmas, N, 1) .* e_mat;

ith_snr = 1;
tic
for it = 1 : num_iter

    e_n  = e_mat(:,it);
    yy_n = y_n + e_n;

    % Store the current SNR
    P_e = (e_n'*e_n) / length(e_n);
    SNRm_vec(it) = 10 * log10( P_y / P_e );

    % Compute the moments of the signal
    s_m = c_m_n * yy_n;

    % Locate Diracs using matrix pencil method
    uu_k          = pencil(s_m, K);
    pp_k          = mod(angle(uu_k), 2*pi);
    ttk_pen(it,:) = sort(T * pp_k / lambda);

    % Locate Diracs using Cadzow+Prony
    s_m = cadzow(s_m, K, cad_it);
    S   = toeplitz(s_m(K+1:end), s_m(K+1:-1:1));
    [~,~,V]       = svd(S);
    h             = V(:,end);
    uu_k          = roots(h);
    pp_k          = mod(angle(uu_k), 2*pi);
    ttk_cad(it,:) = sort(T * pp_k / lambda);
end
toc

% Compute the Cramer-Rao bound
crb_yn = get_crb_yn(phi, t_phi, T, n_vec, t_k(1:K), a_k(1:K), SNRs, P_y);
% crb_sm = get_crb_sm_exp(alpha_vec, T, t_k(1:K), a_k(1:K), c_m_n, SNRs, P_y);

% Compute the actual standard deviation for t_k
std_tk_pen = zeros(length(SNRs), K);
std_tk_cad = zeros(length(SNRs), K);
for it = 1 : length(SNRs)
    ttk_cur          = ttk_pen((it-1)*it_per_snr+1:it*it_per_snr,:);
    std_tk_pen(it,:) = sqrt(sum((ttk_cur - repmat(t_k(1:K), it_per_snr, 1)).^2) / it_per_snr);
    
    ttk_cur          = ttk_cad((it-1)*it_per_snr+1:it*it_per_snr,:);
    std_tk_cad(it,:) = sqrt(sum((ttk_cur - repmat(t_k(1:K), it_per_snr, 1)).^2) / it_per_snr);
end

SNRm_vec1 = repmat(SNRm_vec, K, 1);
ttk_vec1  = ttk_pen(:);

% save('results.mat', 'K', 'SNRs', 'crb_yn', 'std_tk_pen', 'std_tk_cad', 'SNR_min', 'SNR_max')

%%
figure
set(gcf, 'Position', [50 50 420 315])
scatter(SNRm_vec1, ttk_vec1, 25, [0.5 0.5 0.5], 'filled')
hold on
for k = 1 : K
    plot(SNRs, t_k(k)*ones(1,length(SNRs)), 'k', 'LineWidth', 2)
end
axis([min(SNRs) max(SNRs) 0 t_int-0.9*T])
xlabel('$$SNR$$ [dB]', 'FontSize', label_size, 'Interpreter', 'Latex')
ylabel('$$t_k$$', 'FontSize', label_size, 'Interpreter', 'Latex')
% hdl = title(['K = ' num2str(K) ' Diracs, P = ' num2str(P) ', N = ' num2str(N)]);
% set(hdl, 'FontSize', label_size)
grid on
set(gca, 'box', 'on')
set(gca, 'FontSize', axis_size)

%%

load('results.mat', 'K', 'SNRs', 'crb_yn', 'std_tk_pen', 'std_tk_cad', 'SNR_min', 'SNR_max')

label_size = 28;
axis_size  = 24;
text_size  = 28;

for k = 1 : K
    figure
    set(gcf, 'Position', [50+k*250 50 420 315])
    semilogy(SNRs, sqrt(crb_yn(:,k)), '-k', 'LineWidth', 2)
    hold on
    semilogy(SNRs, std_tk_pen(:,k), '--k', 'LineWidth', 2)
    semilogy(SNRs, std_tk_cad(:,k), 'xk', 'LineWidth', 2, 'MarkerSize', 8)
    set(gca, 'XTick', SNR_min:5:SNR_max);
    set(gca, 'YTick', [1e-5,1e-4,1e-3,1e-2,1e-1,1e0]);
    axis([SNR_min SNR_max 5e-4 5e-1])
    grid on
    xlabel('$$SNR$$ [dB]', 'FontSize', label_size, 'Interpreter', 'Latex')
    ylabel(['$$\Delta t_' num2str(k) '$$'], 'FontSize', label_size, 'Interpreter', 'Latex')
    hdl = legend('CRB', ...
                 'Pencil', 'Cadzow');
    set(hdl, 'FontSize', axis_size)
    set(gca, 'FontSize', axis_size)
end

