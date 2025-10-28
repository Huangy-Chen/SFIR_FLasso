function [index1,index2,s] = screeningmax(X,y,p,lambda1,lambda2)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
% opt.n, opt.p, opt.lambda1, opt.lambdatilde,opt.lambda2,
o = y/2;
r = norm(o);
for j = 1: p
xy(j)=norm(X(:,j))*r+abs(X(:,j)'*o);
end
k = 1;
l = 1;
index1 = [];
index2 = [];

for j = 1 : p
    %% active feature selection rule
    if j==1&&(xy(j)>2*(lambda1-lambda2)||xy(j)==2*(lambda1-lambda2))
        index1(k)=1;
        k = k+1;
    else
        if j>1&&j<p&&(xy(j)>2*(lambda1-2*lambda2)||xy(j)==2*(lambda1-2*lambda2))
        index1(k)=j;
        k=k+1;
        else
            if xy(j)>2*(lambda1-lambda2)||xy(j)==2*(lambda1-lambda2)
              index1(k)=p;
            end
        end
    end
    %% successively same coefficient 
    add = sum(xy(1:j)); 
    if add+j*lambda1>lambda2||add+j*lambda1==lambda2
        index2(l) = j;
        l = l+1;
    end
end
s = size(index1,2);
end

