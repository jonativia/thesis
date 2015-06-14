function [phi, t] = generate_e_spline(alpha_vec, T_s, T, mode)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2011
%
% Date        : 21/11/2011
% Supervisor  : Dr Pier Luigi Dragotti
% Author      : Jon Onativia
%
% File        : generate_e_spline.m
% -------------------------------------------------------------------------
% Generate the exponential spline of order P+1 corresponding to a vector of
% alpha values and with a given temporal resolution. The resulting spline 
% is obtained in time domain computing the P convolutions of the P+1 zero 
% order E-splines:
%   phi_a_vec(t) = phi_a_0(t) * phi_a_1(t) * ... * phi_a_N(t)
%
% USAGE:
%  [phi, t] = generate_e_spline(alpha_vec, T_s[, T, mode])
%
% INPUT:
%  - alpha_vec : Vector of P+1 alpha values of the E=spline.
%  - T_s       : Time resolution of the spline.
%  - T         : Optional argument. Scale factor. Default T = 1.
%  - mode      : Optional argument. 'causal', 'symmetric' or 'anticausal'. 
%                Default mode = 'causal'.
%
% OUTPUT:
%  - phi       : Vector of size (P+1)/T + 1 with the values of the
%                E-spline.
%  - t         : Time stamps of the corresponding values of the phi vector.
%

if nargin < 2 || nargin > 4
    error('generate_e_spline:err_arg', 'The number of input arguments is incorrect.')
elseif nargin < 4
    mode = 'causal';
    if nargin < 3
        T = 1;
    end
end

P = length(alpha_vec) - 1;

% Convert alpha_vec into a row vector
alpha_vec = alpha_vec(:).';

% Apply scaling factor
len = (P+1) * T / T_s + 1;
N   = 2^nextpow2(len);
w   = 2*pi/(N*T_s) * (-(N/2) : (N/2 - 1))';

% Build the B-spline in the frequency domain
[X, Y]           = meshgrid(alpha_vec, 1j*w*T);
num              = 1 - exp(X - Y);
denum            = Y - X;
indet_idx        = (num == 0) & (denum == 0);
num(indet_idx)   = 1;
denum(indet_idx) = 1;
phi_w = prod(num./denum, 2);

% Compute the inverse Fourier transform
phi = (T/T_s) * real( ifft( [phi_w(end/2+1:end); phi_w(1:end/2)] ) );
phi = phi(1:len+1);
t  = (0 : len)' * T_s;

if strcmp(mode, 'symmetric')
    t_mid      = (t(end) - t(1)) / 2;
    t          = t - t_mid;
    [~, i_max] = max(phi);
    if phi(i_max) ~= phi(t == 0)
        t = t - t(i_max);
    end
elseif strcmp(mode, 'anticausal')
    phi = phi(end:-1:1);
    t  = -t(end:-1:1);
end
