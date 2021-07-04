clear all;


% Load the matrix P3_gt containing the ground thruth data:
% P3_gt([t t+N t+2*N],:) contains the 3D coordinates of the P points at time t
% (N is the number of frames, P is the number of points)
load('face.mat');
[N, P] = size(P3_gt); N = N/3;

% 2D motion resulting from orthographic projection
% (input to the non-rigid sfm algorithm)
p2_obs = P3_gt(1:2*N, :);

figure
hold on
plot3(P3_gt(1,:), P3_gt(N+1,:), P3_gt(2*N+1,:), 'go')
axis equal


% Choose K, the size of the shape basis
K = 2;
K = 3 * K;

% Form the measurement matrix
X = p2_obs;
A = X(1:N,:);
B = X(N+1:2*N,:);

W = reshape([A(:), B(:)]',2*N, []);

% Centralize
wm = mean(W,2);
W = W - wm*ones(1,size(W,2));

% Decompose the measurement matrix
[U, D, V] = svd(W);

U = U(:, 1:K);
D = D(1:K, 1:K);
V = V(:, 1:K);

Q = U * sqrt(D);
B = sqrt(D) * V';


% Try frame 1
q = Q(1:2,:);
qb = [q(1,1), q(1,2), q(1,3), q(2,1), q(2,2), q(2,3);
        q(1,4), q(1,5), q(1,6), q(2,4), q(2,5), q(2,6)];
    

% Determine the coefficients
[u, d, v] = svd(qb);

% Get the first principal vector
u = u(:, 1);
d = d(1,1);
v = v(1,:);

% Get the weights and the rotation vectors
w = u * sqrt(d);
r = sqrt(d) * v';


B1 = B(1:3,:);
B2 = B(4:6,:);
est = w(1) * B1 + w(2) * B(2);


% Align reconstruction with ground truth
[~, est, ~] = absor(est, [P3_gt(1,:); P3_gt(N+1,:); P3_gt(2*N+1,:)], 'doScale', true);

% Draw both reconstruction and groudtruth on the same figure
figure
hold on
plot3(est(1,:), est(2,:), est(3,:), 'go');
plot3(P3_gt(1,:), P3_gt(N+1,:), P3_gt(2*N+1,:), 'ro');
axis equal

