% A simple example demonstrating Dual Principal Component Pursuit (DPCP).
% The data consist of inliers lying in a linear subspace of high relative 
% dimension, mixed with outliers. The DPCP problem is then solved via
% an IRLS scheme (DPCP_IRLS) and an estimate for the inlier subspace
% produced. Finally, the distance of each point to the subspace is computed.
% This distance can be used for inlier/outlier separation.

% For more details, please see M.C. Tsakiris and R. Vidal
% "Dual Principal Component Pursuit", International Conference on
% Computer Vision", 2015. 

% Copyright @ Manolis C. Tsakiris, 2016

clear all
close all

%% parameters
% ambient dimension
D = 100;
% subspace dimension
dim = 90;
c = D - dim;
% number of inliers
N = 300;
% number of outliers
M = 300;
% noise standard deviation
noise = 0.00;
% DPCP_IRLS parameters
T = 100;
epsilon_J = 10^(-6);
delta = 10^(-9);

% create subspace
V = randn(D);
[Q, R] = qr(V);
span = Q(:,1:dim);
normals = Q(:,dim+1:D);

% sample the subspace
X = span*randn(dim,N)+noise*normals*randn(D-dim,N);
% add outliers
O = randn(D,M);
X_tilde = [X O];
S=X_tilde(:,3);
% DPCP-IRLS
[P_DPCP_IRLS, B_DPCP_IRLS] = DPCP_IRLS(X_tilde,c,delta,T,epsilon_J);

plot(P_DPCP_IRLS,'r');


%DPCP-IRLS
N=[];
for i = 1:600
    if (P_DPCP_IRLS(i)<0.01)
        row=[1,0,0];
        N=[N;row];
    else
        row=[1,1,0];
        N=[N;row];
    end
end

reflectivity = a(4:4:200000);
point_cloud = [x,y,z,reflectivity];
data = pointCloud([x y z],Color=N);
figure(1);
pcshow(data);hold on;

