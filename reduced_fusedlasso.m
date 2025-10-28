function [x1] =  reduced_fusedlasso(A,y,p,index1,s1,index2,s2,lambda1,lambda2, opts);
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
x1 = zeros(p,1);
if s1==0
    x1 = zeros(p,1);
elseif s1==1
        a = A(:,index1)'*y;
        b = norm(A(:,index1))^2;
        x1(index1)=sign(a)*max(abs(a)-lambda1,0)/b;
else
        [x1(index1),~,~] = fusedLeastR(A(:,index1), y, lambda1,lambda2, opts);
end
end

