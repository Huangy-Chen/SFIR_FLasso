clear, clc;

% To evaluate our proposed safe screening rule for the fused Lasso, we employ the fusedLeastR function 
% from the SLEP package (Liu et al., 2009) as a benchmark solver. 

%% Problem
%
%  min  1/2 || A x - y||^2 + lambda1 * ||x||_1 + lambda2 * sum_i |x_i-x_{i+1}|
%
%% SLEP: 
%
% [1] J. Liu, L. Yuan, and J. Ye, An Efficient Algorithm for a Class of
%     Fused Lasso Problems, KDD, 2010.
% [2] J. Liu, S. Ji, and J. Ye, SLEP: Sparse Learning with Efficient Projections. 
%     Arizona State University, 2009.
%
% Download: http://yelabs.net/software/SLEP/


root=cd;
addpath(genpath([root '/SLEP'])); % add the functions in the folder SLEP to the path
                   
randNum=1;
n = 100;
p = 12000;
beta = zeros(p,1);
beta(1) = 2;
beta(3) = 1.5;
beta(5) = 0.8;
beta(8) = 1;
beta(10) = 1.75;
beta(13) = 0.75;
beta(16:50) = 0.3; 
mu = zeros(p,1);
mu(3:7) = 10;
mu(70:90) = 5;
a =round(p/2);
b =round((2*p)/3);
mu(a:b) = -2;
Sigma = speye(p);
% Sigma = zeros(p,p);
% for i =1 : p
%     for j =1 : p
%         inn = abs(i-j);
%         Sigma(i,j) = 0.5^inn;
%     end
% end
A = mvnrnd(mu,Sigma,n);
error = 0.1*rand(n,1);
y = A*beta+error;
lambda2 = 1e-4; 
if norm(A'*y,'inf')==A(:,1)'*y|| norm(A'*y,'inf')==A(:,p)'*y
    lambda1max = lambda2+norm(A'*y,'inf');
else
    lambda1max = 2*lambda2+norm(A'*y,'inf');
end
lambda1 = [0.01:0.01:1].*lambda1max;          % the regularization parameter
%----------------------- Set optional items ------------------------
opts=[];

% Starting point
opts.init=2;        % starting from a zero point

% termination criterion
opts.tFlag=5;       % run .maxIter iterations
opts.maxIter=500;   % maximum number of iterations

% normalization
opts.nFlag=0;       % without normalization

% regularization
opts.rFlag=1;       % the input parameter 'rho' is a ratio in (0, 1)
%opts.rsL2=0.01;    % the squared two norm term

% fused penalty
%opts.fusedPenalty=0.0001;

% line search
opts.lFlag=0;
x1 = cell(100,1);
x2 = cell(100,1);
%----------------------- Run the code LeastR  without screening rule-----------------------
t1=clock;
sparsity = zeros(100,1);
for i = 1: 1: 100
[x1{i}, funVal1, ValueL1]= fusedLeastR(A, y, lambda1(i),lambda2, opts);
sparsity(i) = sum(x1{i}<1e-3);
end
t2=clock;
time1 = etime(t2,t1);

%----------------------- Run the code LeastR with screening rule -----------------------
t3=clock;
x2{100}=zeros(p,1);
s1 = zeros(100,1);
s1(100)=0;
for k = 100 :-1: 2
[index1, s1(k-1),index2,s2(k-1)] = screening(A,y,x2{k},p,lambda1(k-1),lambda1(k),lambda2);
[x2{k-1}] = reduced_fusedlasso(A,y,p,index1,s1(k-1),index2,s2(k-1),lambda1(k-1),lambda2, opts);
end
t4=clock;
time2 = etime(t4,t3);

%----------------------- Run the screening rule -----------------------
t5=clock;
for k = 100 :-1: 2
[index1, s1(k-1),index2,s2(k-1)] = screening(A,y,x2{k},p,lambda1(k-1),lambda1(k),lambda2);
end
t6=clock;
time3 = etime(t6,t5);

speedup = time1/time2;
fprintf('computational time of solver = %f\n',time1);
fprintf('computational time of solver+screening rule = %f\n',time2);
fprintf('computational time of screening rule = %f\n',time3);
fprintf('speedup of the screening rule = %f\n',speedup);
rejectionratio = min((p-s1)./sparsity,1)';
%plot(rejectionratio)
x = [0.01:0.01:1];
plot(x,rejectionratio)