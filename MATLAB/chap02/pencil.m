function u_k = pencil(y, K)

estimate_k = false;

if nargin < 2
    estimate_k = true;
end

% Params
sv_thresh = 0.3;
N         = length(y);

% Construct a square (or almost square) Toeplitz matrix
mid   = ceil(N/2);
Y_toe = toeplitz(y(mid:end), y(mid:-1:1));

% Estimate K
if estimate_k
    sing_vals = svd(Y_toe); 
    sing_vals = sing_vals / sing_vals(1);
    K         = sum(sing_vals > sv_thresh);
end

% Apply the matrix pencil method to retrieve the locations
[U,~,~] = svd(Y_toe);
U       = U(:,1:K);
Z       = pinv(U(1:end-1,:)) * U(2:end,:);
u_k     = eig(Z);

end
