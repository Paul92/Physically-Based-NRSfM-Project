% inputs: 
% bbs: the struct data types created from the bbs code by Brunet 
% ctrlpts (2 X N): BBS control points
% Q (M X 2): feature points on the template
% Qp (M X 2): feature points on the image
% Sw (L X 2): a grid of points where the shwarzian is imposed, usually a 50X50
% er: Smoothing weight

function ctrlpts = optimPanalSchwarz(bbs,ctrlpts,Q,Qp,SQ,er)

% Qp = Qp(:,1:2);
ctrlpts = ctrlpts';
m = size(Q,1);
l = size(ctrlpts,1);

options = optimset('Display','none','TolFun',1e-6,'MaxIter',4,'Jacobian','on');
cf = @(ctrlpts)L_cf(bbs,ctrlpts,Q,Qp,SQ,er,l);
[ctrlpts,RMSR2] = lsqnonlin(cf,ctrlpts(:),[],[],options);

ctrlpts = reshape(ctrlpts,l,2);
RMSR = sqrt(RMSR2/m);
        
end
 
% % cost function for the nonlinear refinement
% function [y, jMat] = L_cf(L,C,Q,Qp,SQ,lambda,l,er)
function [y,jac] = L_cf(bbs,ctrlpts,Q,Qp,SQ,er,l)

m = size(Q,1);
ctrlpts = reshape(ctrlpts,l,2);
coloc=bbs_coloc(bbs, Q(:,1), Q(:,2));
Qt = bbs_eval(bbs, ctrlpts', Q(:,1), Q(:,2))';
yd = (Qp-Qt);
[I J M N jI jJ jM jN] = schwarzian(bbs,ctrlpts',SQ(:,1),SQ(:,2),'noden');

% [JI JJ JM JN] = schwarzian(SQ, L, C, lambda, 'jacobian_noden',MM001,MM10,MM01,MM20,MM02,MM11);


% jMat = [-MM00;er.*JI;er.*JJ;er.*JM;er.*JN];
% jMat = [-MM00;er.*JM;er.*JN];

ys = [I;J;M;N];
% yb = 5*B;
% yb = yb*L;
% ys = [M;N];
ys =er.*ys;
% ys = [yb(:);ys(:)];

y = [yd(:) ; ys(:)];
% disp(sprintf('%f',y'*y));
% fprintf('Iteration, error E:%f S:%f\n', sum(sum(y.^2)),sum(ys.^2));
jacd=[-coloc,zeros(size(coloc));zeros(size(coloc)),-coloc];
jacs=er.*[jI;jJ;jM;jN];
jac=[jacd;jacs];
disp(sprintf('%f\n',y'*y))
end
