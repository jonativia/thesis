function [bP, t] = generate_b_spline(P, T_s, T, mode)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2011
%
% Date        : 21/11/2011
% Supervisor  : Dr Pier Luigi Dragotti
% Author      : Jon Onativia
%
% File        : generate_b_spline.m
% -------------------------------------------------------------------------
% Generate the B-spline of order P.
%   phi(t) = beta_0(t) * beta_0(t) * ... * beta_0(t)
%
% USAGE:
%  [phi, t] = generate_b_spline(P, T_s[, T, mode])
%
% INPUT:
%  - P      : Order of the B-spline.
%  - T_s    : Time resolution of the spline.
%  - T      : Optional argument. Scale factor. Default T = 1.
%  - mode   : Optional argument. 'causal', 'symmetric' or 'anticausal'. 
%             Default mode = 'causal'.
%
% OUTPUT:
%  - phi    : B-Spline.
%  - t      : Time stamps of the corresponding values of the B-spline.
%

if nargin < 2 || nargin > 4
    error('generate_b_spline:err_arg', 'The number of input arguments is incorrect.')
elseif nargin < 4
    mode = 'causal';
    if nargin < 3
        T = 1;
    end
end

% Apply scaling factor
len = (P+1) * T / T_s + 1;
N   = 2^nextpow2(len);
w   = 2*pi/(N*T_s) * (-(N/2) : (N/2 - 1))';

% Build the B-spline in the frequency domain
b0_w = (1 - exp(-1i * w * T)) ./ (1i * w * T);
b0_w(w==0) = 1;
bP_w = b0_w.^(P+1);

% Build the B-spline in the time domain
bP = (T/T_s) * real(ifft([bP_w(end/2+1:end); bP_w(1:end/2)]));
bP = bP(1:len+1);
t  = (0 : len)' * T_s;

% Adjust the temporal shift
if strcmp(mode, 'symmetric')
    [~, i_max] = max(bP);
    t = t - t(i_max);
elseif strcmp(mode, 'anticausal')
    bP = bP(end:-1:1);
    t   = -t(end:-1:1);
end

end
