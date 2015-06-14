clear all
close all
clc

% Font sizes
label_size = 28;
axis_size  = 24;
text_size  = 28;

f_s = 1;
T   = 1/f_s;
fac = 64;
res = T/fac;
t   = (-T:res:6*T)';
f_1 = f_s/4;
f_2 = f_s/7;
x   = 0.8 + .5*cos(2*pi*f_1*t) + .3*cos(2*pi*f_2*t);

% Sampling process
t_n = t(1:fac:end);
x_n = x(1:fac:end);
n   = round(t_n / T);

max_x = max(x);
min_x = min(x);
dyn_r = max_x - min_x;

figure
set(gcf, 'Position', [50 50 420 315])
plot(t, x, 'k', 'LineWidth', 2)
hold on
stem(t_n, x_n, '.k', 'MarkerSize', 40)
xlabel('$$t$$ [s]', 'Interpreter', 'Latex', 'FontSize', text_size)
set(gca, 'FontSize', axis_size)
axis([t(1)+.5*T t(end)-.5*T min_x-.4*dyn_r max_x+.4*dyn_r])

t_s = (-7*T:res:7*T)';
f   = sinc(t_s/T);

max_f = max(f);
min_f = min(f);
dyn_r = max_f - min_f;

figure
set(gcf, 'Position', [250 50 420 315])
plot(t_s, f, 'k', 'LineWidth', 2)
xlabel('$$t$$ [s]', 'Interpreter', 'Latex', 'FontSize', text_size)
set(gca, 'FontSize', axis_size)
axis([t_s(1) t_s(end) min_f-.2*dyn_r max_f+.2*dyn_r])

figure
set(gcf, 'Position', [450 50 420 315])
plot(t, x, '--k', 'LineWidth', 4)
hold on
stem(t_n, x_n, '.k', 'MarkerSize', 40)
for i = 1 : length(n)
    plot(t_s+n(i)*T, x_n(i)*f, 'k', 'LineWidth', 1)
end
xlabel('$$t$$ [s]', 'Interpreter', 'Latex', 'FontSize', text_size)
set(gca, 'FontSize', axis_size)
axis([t(1)+.5*T t(end)-.5*T min_x-.4*dyn_r max_x+.4*dyn_r])


