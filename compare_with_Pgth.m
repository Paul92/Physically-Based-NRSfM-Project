function [P2,err_p] = compare_with_Pgth(P,u,v,q,Pg)
er = 1e-4;
t= 1e-3;
nC = 40;
for i = 1: size(u,1)
    
    q1 = q{1,i};
    umin=min(u(i,:))-0.1;umax=max(u(i,:))+0.1;
    vmin=min(v(i,:))-0.1;vmax=max(v(i,:))+0.1;
    
    bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 3);
    coloc = bbs_coloc(bbs, u(i,:), v(i,:));
    lambdas = er*ones(nC-3, nC-3);
    bending = bbs_bending(bbs, lambdas);
    cpts = (coloc'*coloc + bending) \ (coloc'*P(3*(i-1)+1:3*(i-1)+3,:)');
    ctrlpts = cpts';
    qw = bbs_eval(bbs, ctrlpts, q1(1,:)',q1(2,:)',0,0);
    
    %[qw,~,~] = RegisterToGTH(qw,Pg(3*(i-1)+1:3*(i-1)+3,:));
    [~,qw,~]= absor(qw,Pg(3*(i-1)+1:3*(i-1)+3,:),'doScale',true);
    P2(3*(i-1)+1:3*(i-1)+3,:) = qw;
   
    scale = 1;%max(max(Pg(3*(i-1)+1:3*(i-1)+3,:)')-min(Pg(3*(i-1)+1:3*(i-1)+3,:)'));
    err_p(i,:) = sqrt(mean((Pg(3*(i-1)+1:3*(i-1)+3,:)-P2(3*(i-1)+1:3*(i-1)+3,:)).^2))/scale;
   
end