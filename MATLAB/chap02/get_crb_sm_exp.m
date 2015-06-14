function [crlb] = get_crb_sm_exp(alpha_m, T, t_k, a_k, C, param1, param2)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2011
%
% Date        : 10/01/2012
% Supervisor  : Dr Pier Luigi Dragotti
% Author      : Jon Onativia and Jose Antonio Uriguen
%
% File        : get_crb_sn.m
% -------------------------------------------------------------------------
%  Compute the Cramer-Rao bound for a set of signal-to-noise ratios.
%

if nargin == 7
    SNR = param1(:);
    P_y = param2;
elseif nargin == 6
    sigma2 = param1;
else
    error('get_crb_sm:err_arg', 'The number of input arguments is incorrect.')
end

% Omega in column vector form, parameters in row vector form
alpha_m = alpha_m(:);
t_k = t_k(:).';
a_k = a_k(:).';

% Compute de Fisher information matrix (except the sigma term)
phase = 1/T * alpha_m * t_k;
coef  = 1/T * alpha_m * a_k;
phi   = [coef.*exp(phase) exp(phase)];
R     = C * C';

% Compute the CRLB for each element of the diagonal
s2_Fish_inv = pinv( phi' * pinv(R) * phi );
if nargin == 7
    sigma2 = P_y * 10.^(-SNR/10);
end
crlb = sigma2 * abs(diag(s2_Fish_inv).');
