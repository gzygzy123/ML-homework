function [distances, B] = DPCP_IRLS(X_tilde,c,delta,T,epsilon_J)

% solves min_B ||X^T B||_1 s.t. B^T B=I
% INPUT
% INPUT
% X_tilde  : DxN data matrix of N data points of dimension D
% c        : Dimension of the orthogonal complement of the subspace.
%            To fit a hyperplane use c=1.
% delta    : Avoids division by zero. Typically is set to 10^(-9).
% T        : Maximal number of iterations. Typically set to 100.
% epsilon_J: Convergence accuracy: Typically set to `10^(-6).
% OUTPUT
% distances: Distance of each point to the estimated subspace.
% B        : Dxc matrix containing in its columns an orthonormal basis for
%            the orthogonal complement of the subspace.

% COPYRIGHT @ Manolis C. Tsakiris, 2016

[D, N] = size(X_tilde);
Delta_J = Inf;
k = 0;
w = ones(N,1);
J_old = zeros(N,1);
J_new = zeros(N,1);
while (Delta_J>epsilon_J) && (k<T)
    R_X = X_tilde * diag(w) * X_tilde'; 
    [U,S,V] = svd(R_X);  % U=V  UXVT=SVD() [U,S,V]特指V; 
    B = U(:,D-c+1:D); 
    for j = 1 : N  %更新w(j)
       w(j) = 1/max(delta,norm(B'*X_tilde(:,j)));
       J_new(j) = norm(B'*X_tilde(:,j));
    end
    k = k + 1;
    Delta_J = abs(sum(J_old)-sum(J_new))/(sum(J_old)+10^(-9));
    J_old = J_new;
end
distances = zeros(1,N);
for j = 1 : N
   distances(j) = norm(B' * X_tilde(:,j)); 
end
end

