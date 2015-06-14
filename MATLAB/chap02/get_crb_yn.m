function [crlb] = get_crb_yn(phi, t_phi, T, n_vec, t_k, a_k, param1, param2)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2011
%
% Date        : 10/01/2012
% Supervisor  : Dr Pier Luigi Dragotti
% Author      : Jon Onativia and Jose Antonio Uriguen
%
% File        : get_crb_yn.m
% -------------------------------------------------------------------------
%  Compute the Cramer-Rao bound for a set of signal-to-noise ratios.
%

if nargin == 8
    SNR = param1(:);
    P_y = param2;
elseif nargin == 7
    sigma2 = param1;
else
    error('get_crb_yn:err_arg', 'The number of input arguments is incorrect.')
end

% Signals in column vector form, parameters in row vector form
phi   = phi(:);
t_phi = t_phi(:);
n_vec = n_vec(:);
t_k   = t_k(:).';
a_k   = a_k(:).';

% Sampling resolution, number of parameters and kernel length
T_s = (t_phi(2) - t_phi(1)) / T;
K   = length(t_k);
L   = length(phi);

% Compute the derivative of phi(t) in time
% dphi = [diff(phi); 0] / T_s;

% Compute the derivative of phi(t) in frequency
N    = 2^nextpow2(L);
w    = (2 * pi * [0:N/2-1 -N/2:1:-1]' / (N * T_s));
dphi = ifft(1j * w .* fft(phi, N));
dphi = dphi(1:L);

% Find the indices of the phi vector that correspond to the shifted t_k
N    = length(n_vec);
tt_k = repmat(t_k, N, 1);
nn   = T * repmat(n_vec, 1, K);
tt_k = (tt_k - nn);
idx  = round((tt_k - t_phi(1)) / (T*T_s)) + 1;

% Make the samples that are outside the support point to a zero value
idx(idx<1 | idx>L) = L+1;
phi                = [phi;0];
dphi               = [dphi;0];

% Compute de Fisher information matrix (multiplied by sigma2)
f_coef  = repmat([a_k/T ones(1, K)], N, 1);
grad_f  = f_coef .* [dphi(idx), phi(idx)];
s2_Fish = grad_f.' * grad_f;

% Compute the CRLB for each element of the diagonal
s2_Fish_inv = inv(s2_Fish);
if nargin == 8
    sigma2 = P_y * 10.^(-SNR/10);
end
crlb = sigma2 * abs(diag(s2_Fish_inv).');

