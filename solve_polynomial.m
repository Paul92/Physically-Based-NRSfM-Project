function res = solve_polynomial(eq)
mpol k1
mpol k2
for i = 1:length(eq)  
    pol = [k1^6,k1^5*k2,k1^4*k2^2,k1^3*k2^3,k1^2*k2^4,...
              k1*k2^5,k2^6,k1^5,k1^4*k2,k1^3*k2^2,k1^2*k2^3,k1*k2^4,k2^5,...
              k1^4,k1^3*k2,k1^2*k2^2,k1^1*k2^3,k2^4,k1^3,k1^2*k2,k1^1*k2^2,...
              k2^3,k1^2,k1*k2,k2^2,k1,k2,1]*eq(:,i);
    P = msdp(min(pol));
    [status,obj,M] = msol(P);
    if status ==1
        res(i,:) = double([k1 k2]);
    else
        res(i,:) = [0,0];
    end
end
