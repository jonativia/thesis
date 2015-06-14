clear all
close all

% Font sizes
label_size = 20;
axis_size  = 16;
text_size  = 14;

f_s = 1;
T   = 1/f_s;
fac = 64;
res = T/fac;
R   = 6;
tau = (R+1) * T;

t_s = (-10*T:res:10*T)';
f   = sinc(t_s/T);

max_f = max(f);
min_f = min(f);
dyn_r = max_f - min_f;


figure
set(gcf, 'Position', [50 50 420 315])
plot(t_s-(R/2)*T, f, 'k', 'LineWidth', 2)
hold on
plot(t_s-(R/2-1)*T, f, 'k', 'LineWidth', 2)
plot(t_s+(R/2)*T, f, 'k', 'LineWidth', 2)
plot([-T/3 0 T/3], [.5 .5 .5], '.k', 'MarkerSize', 15)
stem([-tau/2 tau/2], [1.4+.1*T 1.4+.1*T], 'k', 'Marker', 'None', 'LineWidth', 1)
plot([-tau/2-.3*T tau/2+.3*T], [1.4 1.4], 'k', 'LineWidth', 1, 'MarkerSize', 15)
text(-.1*T, 1.55, '$$\tau$$', 'FontSize', text_size, 'Interpreter', 'Latex')
text(-(R/2)*T-.4*T, 1.1, '$$r=0$$', 'FontSize', text_size, 'Interpreter', 'Latex')
text(-(R/2-1)*T-.3*T, 1.1, '$$r=1$$', 'FontSize', text_size, 'Interpreter', 'Latex')
text(+(R/2)*T-.6*T, 1.1, '$$r=R$$', 'FontSize', text_size, 'Interpreter', 'Latex')
xlabel('$$t$$ [s]', 'Interpreter', 'Latex', 'FontSize', label_size)
set(gca, 'FontSize', axis_size)
axis([-0.7*tau +0.7*tau min_f-.3*dyn_r max_f+.6*dyn_r])

