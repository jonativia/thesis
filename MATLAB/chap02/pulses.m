close all
clear all

% Font sizes
label_size = 30;
axis_size  = 26;
text_size  = 30;

t_x = 0:0.05:25;
x   = zeros(size(t_x));
t_k = [2; 10; 18];
a_k = [1; .5; .8];
for it = 1 : length(t_k)
    x = x + a_k(it) * (2 * exp(-3*(t_x-2-t_k(it)).^2) - exp(-2*(t_x-2.5-t_k(it)).^2));
end

max_x = max(x);
min_x = min(x);
dyn_r = max_x - min_x;

figure
set(gcf, 'Position', [50 50 420 315])
plot(t_x/t_x(end), x, 'k', 'LineWidth', 2)
xlabel('$$t$$ [s]', 'Interpreter', 'Latex', 'FontSize', text_size)
set(gca, 'FontSize', axis_size)
axis([0 1 min_x-.2*dyn_r max_x+.2*dyn_r])


tau = 1;
T   = 1 / 16;
T_s = T / 64;
t   = (0 : T_s : tau)';

% Stream of Diracs
itk = [150 180 350 420];
K   = length(itk);
t_k = t(itk).';
a_k = [1.5 1 2 3];

figure
set(gcf, 'Position', [250 50 420 315])
stem(t_k, a_k, '^k', 'fill', 'LineWidth', 1, 'MarkerSize', 10)
axis([t(1) t(end) -.5 4])
xlabel('$$t$$ [s]', 'Interpreter', 'Latex', 'FontSize', text_size)
set(gca, 'FontSize', axis_size)

% Piecewise sinusoidals
amps   = [.9  .5  .7  .85 .42 .75 .81 .32];
frecs  = [12  35  23  18  15  18  19  5];
t_disc = [.08 .24 .35 .49 .68 .75 .89];
t_disc = [-Inf t_disc Inf];
x_sin  = zeros(size(t));
for ith = 1 : length(t_disc)-1
    idx = (t>=t_disc(ith) & t<t_disc(ith+1));
    x_sin(idx) = amps(ith) * sin(2 * pi * frecs(ith) * t(idx)); 
end
figure
set(gcf, 'Position', [450 50 420 315])
plot(t, x_sin, 'k', 'LineWidth', 2)
axis([t(1) t(end) -1.2 1.2])
xlabel('$$t$$ [s]', 'Interpreter', 'Latex', 'FontSize', text_size)
set(gca, 'FontSize', axis_size)

% Stream of decaying exponentials
alpha = 20;
x_dec = zeros(size(t));
for k = 1:length(t_k)
    x_dec(t>=t_k(k)) = x_dec(t>=t_k(k)) + ...
        a_k(k) * exp(-(t(t>=t_k(k)) - t_k(k)) * alpha);
end

figure
set(gcf, 'Position', [650 50 420 315])
plot(t, x_dec, 'k', 'LineWidth', 2)
axis([t(1) t(end) -.5 4])
xlabel('$$t$$ [s]', 'Interpreter', 'Latex', 'FontSize', text_size)
set(gca, 'FontSize', axis_size)
