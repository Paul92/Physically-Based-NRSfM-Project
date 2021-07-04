function [Pgth,Ngth,q_n,I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb] = create_tshirt_dataset(idx,n, P_gth, q, K, par)

for i = 1:n
    P2_n{i} = P_gth{idx(i)}; % Ground truth 3D points
    q_n{i} = q{idx(i)}; % Image points corresponding to ground truth
end

% create Ground truth
er = 1e-5;
t = 1e-3; % domain padding
nC = 40; % number of warp control points
for i=1:n
    umin = min(q_n{i}(1,:))-t; umax = max(q_n{i}(1,:))+t;
    vmin = min(q_n{i}(2,:))-t; vmax = max(q_n{i}(2,:))+t;
    bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 3);
    coloc = bbs_coloc(bbs, q_n{i}(1,:), q_n{i}(2,:));
    lambdas = er*ones(nC-3, nC-3);
    bending = bbs_bending(bbs, lambdas);

    cpts = (coloc'*coloc + bending) \ (coloc'*P2_n{i}');

    ctrlpts = cpts';
    dqu = bbs_eval(bbs, ctrlpts, q_n{i}(1,:)',q_n{i}(2,:)',1,0);
    dqv = bbs_eval(bbs, ctrlpts, q_n{i}(1,:)',q_n{i}(2,:)',0,1);
    nu = [dqu(1,:)./sqrt(sum(dqu.^2));...
        dqu(2,:)./sqrt(sum(dqu.^2));...
        dqu(3,:)./sqrt(sum(dqu.^2))];
    
    nv = [dqv(1,:)./sqrt(sum(dqv.^2));...
        dqv(2,:)./sqrt(sum(dqv.^2));...
        dqv(3,:)./sqrt(sum(dqv.^2))];
    
    nn = -cross(nu,nv);
    Pgth(3*(i-1)+1:3*(i-1)+3,:) = P2_n{i};
    Ngth(3*(i-1)+1:3*(i-1)+3,:) = [nn(1,:)./sqrt(sum(nn.^2));...
        nn(2,:)./sqrt(sum(nn.^2));...
        nn(3,:)./sqrt(sum(nn.^2))];
end

% make a grid of points on the images and compute schwarzian warps
er= 1e-4;
p = 20; %size of grid
I2u = zeros(n-1,p*p);
I2v = zeros(n-1,p*p);
J21a = zeros(n-1,p*p);
J21b = zeros(n-1,p*p);
J21c = zeros(n-1,p*p);
J21d = zeros(n-1,p*p);
J12a = zeros(n-1,p*p);
J12b = zeros(n-1,p*p);
J12c = zeros(n-1,p*p);
J12d = zeros(n-1,p*p);
H21uua = zeros(n-1,p*p);
H21uub = zeros(n-1,p*p);
H21uva = zeros(n-1,p*p);
H21uvb = zeros(n-1,p*p);
H21vva = zeros(n-1,p*p);
H21vvb = zeros(n-1,p*p);

% make a grid of points
nC = 40;
q1 = q_n{1}(1:2,:);  
umin = min(q1(1,:))-t; umax = max(q1(1,:))+t;
vmin = min(q1(2,:))-t; vmax = max(q1(2,:))+t;
bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 2);
coloc = bbs_coloc(bbs, q1(1,:), q1(2,:));
lambdas = er*ones(nC-3, nC-3);
bending = bbs_bending(bbs, lambdas);
[xv,yv]=meshgrid(linspace(bbs.umin,bbs.umax,p),linspace(bbs.vmin,bbs.vmax,p));
I1u = xv(:)';
I1v = yv(:)';

for i = 2:n
    q2 = q_n{i}(1:2,:);
    cpts = (coloc'*coloc + bending) \ (coloc'*q2');

    ctrlpts = cpts';
    q1 =  bbs_eval(bbs,ctrlpts,q_n{1}(1,:)',q_n{1}(2,:)',0,0);
    error=sqrt(mean((q1(1,:)-q_n{i}(1,:)).^2+(q1(2,:)-q_n{i}(2,:)).^2));

    %disp(fprintf('[ETA] Internal Rep error = %f',error));
    q =  bbs_eval(bbs,ctrlpts,I1u',I1v',0,0);
    I2u(i-1,:) = q(1,:);
    I2v(i-1,:) = q(2,:);
    
    %Visualize Point Registration Error
%     figure;
%     plot(q_n{i}(1,:),q_n{i}(2,:),'ro');
%     hold on;
%     plot(q(1,:),q(2,:),'b*');
%     %mesh(reshape(q(1,:),size(xv)),reshape(q(2,:),size(xv)),zeros(size(xv)));
%     axis equal
%     hold off;
end

% calculate Eta_21 derivatives using schwarzian warps
for i = 2:n
    q1 = [I1u;I1v];  q2 = [I2u(i-1,:);I2v(i-1,:)];
    umin = min(q2(1,:))-t; umax = max(q2(1,:))+t;
    vmin = min(q2(2,:))-t; vmax = max(q2(2,:))+t;
    bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 2);
    coloc = bbs_coloc(bbs, q2(1,:), q2(2,:));
    lambdas = er*ones(nC-3, nC-3);
    bending = bbs_bending(bbs, lambdas);
    cpts = (coloc'*coloc + bending) \ (coloc'*q1');
    ctrlpts = cpts';
     
    % get control points for j to i warp
    cpts = (coloc'*coloc + bending) \ (coloc'*q1(1:2,:)');
    ctrlpts = cpts';
    [xv,yv]=meshgrid(linspace(umin,umax,p),linspace(vmin,vmax,p));
    
    % COMMENT THESE 2 TO DISABLE SCHWARPS
%       ctrlpts = optimPanalSchwarz(bbs,ctrlpts,q2',q1',[xv(:),yv(:)],par(i));
%       ctrlpts = ctrlpts';
    
    
    dqu = bbs_eval(bbs, ctrlpts, q2(1,:)',q2(2,:)',1,0);
    dqv = bbs_eval(bbs, ctrlpts, q2(1,:)',q2(2,:)',0,1);
    dquv = bbs_eval(bbs,ctrlpts,q2(1,:)',q2(2,:)',1,1);
    dquu = bbs_eval(bbs, ctrlpts, q2(1,:)',q2(2,:)',2,0);
    dqvv = bbs_eval(bbs, ctrlpts, q2(1,:)',q2(2,:)',0,2);
    
    J21a(i-1,:) = dqu(1,:);
    J21b(i-1,:) = dqu(2,:);
    J21c(i-1,:) = dqv(1,:);
    J21d(i-1,:) = dqv(2,:);
    
    J12a(i-1,:) = dqv(2,:)./(dqu(1,:).*dqv(2,:)-dqv(1,:).*dqu(2,:));
    J12b(i-1,:) = -dqu(2,:)./(dqu(1,:).*dqv(2,:)-dqv(1,:).*dqu(2,:));
    J12c(i-1,:) = -dqv(1,:)./(dqu(1,:).*dqv(2,:)-dqv(1,:).*dqu(2,:));
    J12d(i-1,:) = dqu(1,:)./(dqu(1,:).*dqv(2,:)-dqv(1,:).*dqu(2,:));
    
    H21uua(i-1,:) = dquu(1,:);
    H21uub(i-1,:) = dquu(2,:);
    
    H21uva(i-1,:) = dquv(1,:);
    H21uvb(i-1,:) = dquv(2,:);
    
    H21vva(i-1,:) = dqvv(1,:);
    H21vvb(i-1,:) = dqvv(2,:);
    
end

