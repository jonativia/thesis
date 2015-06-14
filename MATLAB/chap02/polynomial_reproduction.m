clear all
close all

% Font sizes
label_size = 30;
axis_size  = 26;
text_size  = 30;

% Sampling period in seconds
T = 1;

% "Continuous" time resolution
T_s = T/256;

% Degree of the B-spline
P = 3;

% Number of samples
N = 13;

% Obtain the sampling kernel
[phi, t_phi] = generate_b_spline(P, T_s, T);

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
c_m_n = get_c_m_n_poly(P, n_vec, phi, t_phi, T);

for ith_m = 1:P+1
    m = ith_m - 1;
    
    % Original polynomial
    f = (t / T).^m;

    % Polynomial reproduction
    f_rec = zeros(length(t), length(n_vec)+1);
    for i = 1:length(n_vec)
        [~, idx_f, idx_phi] = intersect(t, t_phi+n_vec(i)*T);
        f_rec(idx_f,i+1)    = c_m_n(ith_m, i) * phi(idx_phi);

        f_rec(:,1) = f_rec(:,1) + f_rec(:,i+1);
    end

    % Dynamic range of the reproduced polynomial
    dyn_range = max(f_rec(:,1)) - min(f_rec(:,1));
    
    % Measure the error
    phi_support = length(phi);
    idx         = (1+phi_support:length(f)-phi_support)';
    rec_er      = f(idx) - f_rec(idx,1);
    MSE         = (rec_er' * rec_er) / length(rec_er);
    disp(['MSE = ' num2str(MSE)])

    figure
    set(gcf, 'Position', [50+100*ith_m 50 420 315])
    plot(t, f_rec(:,1), 'k', 'LineWidth', 4)
    hold on
    for i = 1:length(n_vec)
        plot(t, f_rec(:,i+1), '-k')
    end
    plot(t, f, '--k', 'LineWidth', 2)
    hdl = xlabel('$$t$$ [s]', 'Interpreter', 'Latex');
    set(hdl, 'FontSize', label_size)
    set(gca, 'FontSize', axis_size)
    axis([t(1) t(end) min(f_rec(:,1))-.3*dyn_range max(f_rec(:,1))+.3*dyn_range])
    
end

