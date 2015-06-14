function c_m_n = get_c_m_n_poly(P, n, phi, t_phi, T, t_0)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2011
%
% Date        : 23/02/2012
% Supervisor  : Dr Pier Luigi Dragotti
% Author      : Jon Onativia
%
% File        : get_c_m_n_poly.m
% -------------------------------------------------------------------------
% Compute the c_m_n coefficients to reproduce polynomials of order up to
% P:
%   t^m = sum_n ( c_m_n * phi(t-n) ), where m = 0, 1, ..., P
%
% USAGE:
%  c_m_n = get_c_m_n_exp(P, n, phi, t[, T, t_0])
%
% INPUT:
%  - alpha_m : Vector of size M with the parameters of the exponentials to 
%              be reproduced.
%  - n       : Vector of size N with the values where the summation will be
%              evaluated.
%  - phi     : Exponential reproducing kernel.
%  - t_phi   : Time stamps of the kernel.
%  - T       : Optional argument. Scale factor. Default T = 1.
%  - t_0     : Optional argument. t value where c_m_0 will be evaluated. Default t_0 = 0.
%
% OUTPUT:
%  - c_m_n   : Coefficients to reproduce the exponentials.
%

if nargin < 4 || nargin > 6
    error('generate_e_spline:err_arg', 'The number of input arguments is incorrect.')
elseif nargin < 6
    t_0 = 1;
    if nargin < 5
        T = 1;
    end
end

% Rearrange the arguments (n row vector, alpha_m column vector)
n       = n(:).';
n_len   = length(n);
T_s     = t_phi(2) - t_phi(1);

% Output matrix
c_m_n = zeros(P+1, n_len);

% Kernel's boundaries
t_1 = t_phi(1) / T;
t_2 = t_phi(end) / T;

% Compute c_m_0 vector
l     = ceil(t_0/T - t_2) : floor(t_0/T - t_1);
idx   = round( (t_0 - T * (t_1 + l)) / T_s) + 1;
phi_l = phi(idx);
sum_p = sum(phi_l);
c_m_0 = zeros(P+1, 1);
for m = 0 : P
    c_m_0(m+1) = (t_0 / T) .^ m;
    for k = 0 : m-1
        n_phi      = nchoosek(m, k) * c_m_0(k+1) * (l.^(m-k) * phi_l);
        c_m_0(m+1) = c_m_0(m+1) - n_phi;
    end
    c_m_0(m+1) = c_m_0(m+1) / sum_p;
end

% From c_m_0 vector, compute the remaining c_m_n coefficients for n ~= 4
n_idx = find(n ~= 0);
for ith_n = n_idx
    cur_n = n(ith_n);
    for m = 0 : P
        for k = 0 : m
            c_m_n(m+1,ith_n) = c_m_n(m+1,ith_n) + nchoosek(m, k) * cur_n^(m-k) * c_m_0(k+1);
        end
    end
end

% Store the c_m_0 values in the c_m_n matrix
n_0          = (n == 0);
c_m_n(:,n_0) = c_m_0;
