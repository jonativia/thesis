function y = cadzow(y, K, iter)

% Construct a square (or almost square) Toeplitz matrix
mid   = ceil(length(y)/2);
Y_toe = toeplitz(y(mid:end), y(mid:-1:1));
[m,n] = size(Y_toe);

for it = 1 : iter
    
    % Force the rank of the matrix to be K
    [U,S,V] = svd(Y_toe);
    Y_toe   = U(:,1:K) * S(1:K,1:K) * V(:,1:K)';
    
    % Average the diagonal elements to obtain a Toeplitz matrix
    c = zeros(m, 1);
    r = zeros(1, n);
    c(1) = mean(diag(Y_toe));
    r(1) = c(1);
    for ith_row = 2 : m
        c(ith_row) = mean(diag(Y_toe,1-ith_row));
    end
    for ith_col = 2 : n
        r(ith_col) = mean(diag(Y_toe,ith_col-1));
    end
    Y_toe = toeplitz(c, r);
end

y = [Y_toe(1,end:-1:1).'; Y_toe(2:end,1)];

end
