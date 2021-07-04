%INPUT
% P : feature points
% C : Control Center
% L : Warp Parameters
% lambda: Inernal Regulrization Parameter

%OUTPUT
% I,J,M,N : The differential invariants.


% function [I J M N] = schwarzian(P, L, C, lambda,choice,MM00,MM10,MM01,MM20,MM02,MM11)
function [I J M N jI jJ jM jN] = schwarzian(bbs,ctrlpts,X2,Y2,choice)



 DAdetabydu1  = bbs_eval(bbs, ctrlpts, X2, Y2, 1,0)';
 DAdetabydv1  = bbs_eval(bbs, ctrlpts, X2, Y2, 0,1)';
 DAdetabydu2  = bbs_eval(bbs, ctrlpts, X2, Y2, 2,0)';
 DAdetabydv2  = bbs_eval(bbs, ctrlpts, X2, Y2, 0,2)';
 DAdetabydudv = bbs_eval(bbs, ctrlpts, X2, Y2, 1,1)';
coloc_du = bbs_coloc_deriv(bbs, X2, Y2,1,0);
coloc_dv = bbs_coloc_deriv(bbs, X2, Y2,0,1);
coloc_duu= bbs_coloc_deriv(bbs, X2, Y2,2,0);
coloc_duv= bbs_coloc_deriv(bbs, X2, Y2,1,1);
coloc_dvv= bbs_coloc_deriv(bbs, X2, Y2,0,2);
 if(strcmp(choice,'den'))
 
      I = (DAdetabydu2(:,1).*DAdetabydu1(:,2) - DAdetabydu2(:,2).*DAdetabydu1(:,1)) ...
                ./ (DAdetabydu1(:,1).*DAdetabydv1(:,2) - DAdetabydu1(:,2).*DAdetabydv1(:,1));

      J = (DAdetabydv2(:,2).*DAdetabydv1(:,1) - DAdetabydv2(:,1).*DAdetabydv1(:,2)) ...
                ./ (DAdetabydu1(:,1).*DAdetabydv1(:,2) - DAdetabydu1(:,2).*DAdetabydv1(:,1));

      M = (DAdetabydu2(:,1).*DAdetabydv1(:,2) - DAdetabydu2(:,2).*DAdetabydv1(:,1) + ...
                2.*(DAdetabydudv(:,1).*DAdetabydu1(:,2) - DAdetabydudv(:,2).*DAdetabydu1(:,1))) ...
                ./ 3.*((DAdetabydu1(:,1).*DAdetabydv1(:,2) - DAdetabydu1(:,2).*DAdetabydv1(:,1)));

      N = (DAdetabydv2(:,2).*DAdetabydu1(:,1) - DAdetabydv2(:,1).*DAdetabydu1(:,2) + ...
                2.*(DAdetabydudv(:,2).*DAdetabydv1(:,1) - DAdetabydudv(:,1).*DAdetabydv1(:,2))) ...
                ./ 3.*((DAdetabydu1(:,1).*DAdetabydv1(:,2) - DAdetabydu1(:,2).*DAdetabydv1(:,1)));

                       
%   elseif(strcmp(choice,'jacobian_noden'))
% %        I = (DAdetabydu2(:,1).*DAdetabydu1(:,2) - DAdetabydu2(:,2).*DAdetabydu1(:,1));
% %        I = (DAdetabydu2(:,1) * MM10 - MM20)
%          coldu2_1 = DAdetabydu2(:,1)*ones(1,size(MM00,2));
%          coldu1_1 = DAdetabydu1(:,1)*ones(1,size(MM00,2));
%          
%          coldu1_2 = DAdetabydu1(:,2)*ones(1,size(MM00,2));
%          coldu2_2 = DAdetabydu2(:,2)*ones(1,size(MM00,2));
%          
%          coldv2_1 = DAdetabydv2(:,1)*ones(1,size(MM00,2));
%          coldv1_1 = DAdetabydv1(:,1)*ones(1,size(MM00,2));
% 
%          coldv1_2 = DAdetabydv1(:,2)*ones(1,size(MM00,2));
%          coldv2_2 = DAdetabydv2(:,2)*ones(1,size(MM00,2));
%          
%          coldudv_1 = DAdetabydudv(:,1)*ones(1,size(MM00,2));
%          coldudv_2 = DAdetabydudv(:,2)*ones(1,size(MM00,2));
%          
%          hf1 = 1:size(MM00,1)/2; 
%          hf2 = size(MM00,1)/2+1:size(MM00,1); 
%          
%          I = coldu1_2.*MM20(hf1,:) + coldu2_1.*MM10(hf2,:) - coldu1_1.*MM20(hf2,:) - coldu2_2.*MM10(hf1,:);
%          
%          J = coldv1_1.*MM02(hf2,:) + coldv2_2.*MM01(hf1,:) - coldv1_2.*MM02(hf1,:) - coldv2_1.*MM01(hf2,:);
%          
%          M = coldv1_2.*MM20(hf1,:) + coldu2_1.*MM01(hf2,:) - coldv1_1.*MM20(hf2,:) - coldu2_2.*MM01(hf1,:) + ...
%                 2.*(coldu1_2.*MM11(hf1,:) + coldudv_1.*MM10(hf2,:) - coldu1_1.*MM11(hf2,:) - coldudv_2.*MM10(hf1,:)); 
% 
%          N = coldu1_1.*MM02(hf2,:) + coldv2_2.*MM10(hf1,:) - coldu1_2.*MM02(hf1,:) - coldv2_1.*MM10(hf2,:) + ...
%                 2.*(coldv1_1.*MM11(hf2,:) + coldudv_2.*MM01(hf1,:) - coldv1_2.*MM11(hf1,:) - coldudv_1.*MM01(hf2,:)); 
% 

  elseif (strcmp(choice, 'noden'))
        
      I = (DAdetabydu2(:,1).*DAdetabydu1(:,2) - DAdetabydu2(:,2).*DAdetabydu1(:,1));

      J = (DAdetabydv2(:,2).*DAdetabydv1(:,1) - DAdetabydv2(:,1).*DAdetabydv1(:,2));

      M = (DAdetabydu2(:,1).*DAdetabydv1(:,2) - DAdetabydu2(:,2).*DAdetabydv1(:,1) + ...
                2.*(DAdetabydudv(:,1).*DAdetabydu1(:,2) - DAdetabydudv(:,2).*DAdetabydu1(:,1)));

      N = (DAdetabydv2(:,2).*DAdetabydu1(:,1) - DAdetabydv2(:,1).*DAdetabydu1(:,2) + ...
                2.*(DAdetabydudv(:,2).*DAdetabydv1(:,1) - DAdetabydudv(:,1).*DAdetabydv1(:,2)));
            
    D12u=DAdetabydu1';
    D12v=DAdetabydv1';
    D12uv=DAdetabydudv';
    D12uu=DAdetabydu2';
    D12vv=DAdetabydv2';
    nparam=size(ctrlpts,2);
    zerosm=zeros(size(coloc_du));
    jI=repmat(D12u(2,:)',1,nparam*2).*[coloc_duu,zerosm]+repmat(D12uu(1,:)',1,nparam*2).*[zerosm,coloc_du];
    jI=jI-(repmat(D12u(1,:)',1,nparam*2).*[zerosm,coloc_duu]+repmat(D12uu(2,:)',1,nparam*2).*[coloc_du,zerosm]);

    jJ=-repmat(D12v(2,:)',1,nparam*2).*[coloc_dvv,zerosm]-repmat(D12vv(1,:)',1,nparam*2).*[zerosm,coloc_dv];
    jJ=jJ+(repmat(D12v(1,:)',1,nparam*2).*[zerosm,coloc_dvv]+repmat(D12vv(2,:)',1,nparam*2).*[coloc_dv,zerosm]);      
    
    jM=repmat(D12v(2,:)',1,nparam*2).*[coloc_duu,zerosm]+repmat(D12uu(1,:)',1,nparam*2).*[zerosm,coloc_dv];
    jM=jM-(repmat(D12v(1,:)',1,nparam*2).*[zerosm,coloc_duu]+repmat(D12uu(2,:)',1,nparam*2).*[coloc_dv,zerosm]);
    jM=jM+2.*(repmat(D12u(2,:)',1,nparam*2).*[coloc_duv,zerosm]+repmat(D12uv(1,:)',1,nparam*2).*[zerosm,coloc_du]);
    jM=jM-2.*(repmat(D12u(1,:)',1,nparam*2).*[zerosm,coloc_duv]+repmat(D12uv(2,:)',1,nparam*2).*[coloc_du,zerosm]);

    jN=-repmat(D12u(2,:)',1,nparam*2).*[coloc_dvv,zerosm]-repmat(D12vv(1,:)',1,nparam*2).*[zerosm,coloc_du];
    jN=jN+(repmat(D12u(1,:)',1,nparam*2).*[zerosm,coloc_dvv]+repmat(D12vv(2,:)',1,nparam*2).*[coloc_du,zerosm]);
    jN=jN-2.*(repmat(D12v(2,:)',1,nparam*2).*[coloc_duv,zerosm]+repmat(D12uv(1,:)',1,nparam*2).*[zerosm,coloc_dv]);
    jN=jN+2.*(repmat(D12v(1,:)',1,nparam*2).*[zerosm,coloc_duv]+repmat(D12uv(2,:)',1,nparam*2).*[coloc_dv,zerosm]);     
    
    
  end
 end

%  
%          coldu2_1 = [DAdetabydu2(:,1) ; zeros(size(DAdetabydu2(:,1)))]*ones(1,size(MM00,2));
%          coldu1_1 = [DAdetabydu1(:,1) ; zeros(size(DAdetabydu1(:,1)))]*ones(1,size(MM00,2));
%          
%          coldu1_2 = [zeros(size(DAdetabydu1(:,1)));DAdetabydu1(:,2)]*ones(1,size(MM00,2));
%          coldu2_2 = [zeros(size(DAdetabydu2(:,1)));DAdetabydu2(:,2)]*ones(1,size(MM00,2));
%          
%          coldv2_1 = [DAdetabydv2(:,1) ; zeros(size(DAdetabydv2(:,1)))]*ones(1,size(MM00,2));
%          coldv1_1 = [DAdetabydv1(:,1) ; zeros(size(DAdetabydv1(:,1)))]*ones(1,size(MM00,2));
% 
%          coldv1_2 = [zeros(size(DAdetabydv1(:,2)));DAdetabydv1(:,2)]*ones(1,size(MM00,2));
%          coldv2_2 = [zeros(size(DAdetabydv2(:,2)));DAdetabydv2(:,2)]*ones(1,size(MM00,2));
%          
%          coldudv_1 = [DAdetabydudv(:,1);zeros(size(DAdetabydudv(:,1)))]*ones(1,size(MM00,2));
%          coldudv_2 = [zeros(size(DAdetabydudv(:,2)));DAdetabydudv(:,2)]*ones(1,size(MM00,2));
