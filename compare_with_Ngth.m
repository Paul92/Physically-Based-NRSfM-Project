function [N,err_n] = compare_with_Ngth(P,q,Ng)
er = 1e-4;
nC = 40;
for i = 1: size(P,1)/3
    
    q1 = q{1,i};
    umin=min(q1(1,:))-0.1;umax=max(q1(1,:))+0.1;
    vmin=min(q1(2,:))-0.1;vmax=max(q1(2,:))+0.1;
    bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 3);
    coloc = bbs_coloc(bbs, q1(1,:), q1(2,:));
    lambdas = er*ones(nC-3, nC-3);
    bending = bbs_bending(bbs, lambdas);
    cpts = (coloc'*coloc + bending) \ (coloc'*P(3*(i-1)+1:3*(i-1)+3,:)');
    ctrlpts = cpts';
    qw = bbs_eval(bbs, ctrlpts, q1(1,:)',q1(2,:)',0,0);
    error = sqrt(mean(sum((qw-P(3*(i-1)+1:3*(i-1)+3,:)).^2)));
    dqu = bbs_eval(bbs, ctrlpts, q1(1,:)',q1(2,:)',1,0);
    dqv = bbs_eval(bbs, ctrlpts, q1(1,:)',q1(2,:)',0,1);
    nu = [dqu(1,:)./sqrt(sum(dqu.^2));...
        dqu(2,:)./sqrt(sum(dqu.^2));...
        dqu(3,:)./sqrt(sum(dqu.^2))];
    
    nv = [dqv(1,:)./sqrt(sum(dqv.^2));...
        dqv(2,:)./sqrt(sum(dqv.^2));...
        dqv(3,:)./sqrt(sum(dqv.^2))];
    
    nn = -cross(nu,nv);
    N(3*(i-1)+1:3*(i-1)+3,:) = [nn(1,:)./sqrt(sum(nn.^2));...
        nn(2,:)./sqrt(sum(nn.^2));...
        nn(3,:)./sqrt(sum(nn.^2))];
    
    err_n(i,:) = acosd(sum(N(3*(i-1)+1:3*(i-1)+3,:).*Ng(3*(i-1)+1:3*(i-1)+3,:)));
    
end