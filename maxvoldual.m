%   function [v1, West, best_theta, iter, Y, C] = maxvoldual(X,r,lambda,options)
% 
%    This code solves simplex-structured matrix factorization (SSMF) 
%    via a dual approach. Given the input matrix, X, and a factorization
%    rank r, it looks for W and H such that WH approximates X and H is column 
%    stochastic. 
%    To do so, it first reduces the dimension of the problem: 
%        Y = C' (X - v e') where e is the vector of all ones, v is in conv(X)
%                                C' contains the first r singular vectors of 
%                                X - v e'. 
%    Then it solves the maximum-volume dual problems: 
%        max_{Z,Theta,Delta}  det(Z)^2 - lambda ||Delta||^2
%                  such that  Z = [Theta; e'] and Y' Theta <= 1 + Delta, 
%    where Theta represents the polar of conv(W) (in the reduced spase) 
%    whose volume is maximized. 
% 
% See the paper "Dual Simplex Volume Maximization for Simplex-Structured 
% Matrix Factorization", by M. Abdolali, G. Barbarino and N. Gillis, 2024.
% 
% ****** Input ******
% X      :  the input matrix
% r      :  the rank of the sought approximation
% lambda :  the regularization parameter
% ---Options---
% .maxiter    : the maximum number of iterations performed by the algorithm
%             -default = 20.
% .epsilon    : the tolerance level for convergence
%             -default = 1e-2.
% .num_workers: number of parallelized solutions used
%             -default = 5.
% .timelimit  : maximum alloted time for the outer loops
%             -default = 300.
% 
% ****** Output ******
% v1          :    estimated center vector
% West        :    estimated W
% best_theta  :    final Theta at the convergence with maximum dual volume
% iter        :    number of iterations till convergence
% Y           :    projected points in r-1 dimensions
% C           :    projection matrix

function [W_est, H_est, Y_res] = maxvoldual(X,r,lambda,options)
tic; 
if nargin <= 3
    options = [];
end
if ~isfield(options,'maxiter')
    options.maxiter = 20;
end
if ~isfield(options,'epsilon')
    options.epsilon = 1e-2;
end
if ~isfield(options,'num_workers')
    options.num_workers = 5;
end
if ~isfield(options,'timelimit')
    options.timelimit = 300;
end
% pre-processing
MAX =  max(max(X));
X = X / MAX ; %scaling for numerical stability
v = mean(X,2);
if issparse(X)
    [m,n]=size(X); %if X is sparse, compute svds implicitly
    [C,~,~]=svds(@afun,[m,n],r-1);
else
    Y = X - v; %otherwise compute svds normally
    [C,~,~]=svds(Y,r-1);
end
% initialization
iter = 0; %initializing iteration counter
v1 = 0; % initializing v_{k-1} (keeps center vector of previous iteration)
nn = 0; % a counter for # of endmembers in SNPA initalization if average initialization goes wrong
theta = cell(options.num_workers,1); % Theta_i for i=1,...num_workers 
ignore = cell(options.num_workers,1); % ignore_i for i=1,...num_workers (ignoring i-th solution or not)
delta = cell(options.num_workers,1); % delta_i for i=1,...num_workers (estimated noise matrix)
flag = cell(options.num_workers,1); % flag_i for i=1,...num_workers (optimization was successful or not)
z = cell(options.num_workers,1); % z_i for i=1,...num_workers
for i = 1: options.num_workers
    z{i}=[randn(r-1,r);ones(1,r)]; %initialize z_i
end
CtX = C'*X;
O = ones(r-1,1);
cro = nchoosek(1:r,r-1); % for calculating intersection later
outeriter = 1;
West = []; 
% main loop
while norm(v-v1,'fro')/norm(v1,'fro') > options.epsilon ...
        && iter < options.maxiter ...
        && toc <= options.timelimit  
   % fprintf('Outer iteration %1.0d started. \n',outeriter); 
    % projection
    v1 = v;
    Y = CtX - C'*v;
    iter = iter + 1;
    % parallel computation of Z_i for i = 1: num_workers
    parfor i=1:options.num_workers
        [z{i},theta{i},delta{i},flag{i}] = algorithm2_update_theta(Y,r,lambda,z{i}); %solves a QP (refer to the paper for more info)
        ignore{i} = ~flag{i} || any(is_correct(normc(theta{i}))<0.01); % check if optimization was successful
    end
    % select the best candidate till now (maximum dual volume)
    best_theta = []; 
    best = -Inf; 
    for i = 1: options.num_workers
        if ~ignore{i} %if optimization for i-th solution was successful
            vol = (det(z{i}))^2 - lambda*sum(delta{i}(:).^2); % evaluate the objective function for i-th candidate
            if vol > best
                best_theta = theta{i}; %update the best candidate
                best = vol;
            end
        else
            z{i}=[randn(r-1,r);ones(1,r)]; % if optimization was unseccessful, change initialization for next iteration
        end
    end
    % if all candidates failed, use alternative initialization
    if isempty(best_theta)
        nn = nn + 1;
        [J,~] = SNPA(X,nn*r); %estimate endmembers with SNPA
        v = mean(X(:,J),2);
        if issparse(X)
            [m,n]=size(X); %if X is sparse, compute svds implicitly
            [C,~,~]=svds(@afun,[m,n],r-1);
        else
            Y = X - v;
            [C,~,~]=svds(Y,r-1); %if not sparse, compute svds normally
        end
        CtX = C'*X; %update CtX for next iterations
        v1 = 0;
    else
        % find intersections (W)
        W_e = [];
        for i=1:size(cro,1) %for each r-1 facets compute intersection
            c = cro(i,:);
            U = best_theta(:,c);
            coef = U' \ O;
            W_e = [W_e coef];
        end
        % update mean vector
        West = C * W_e + v;
        v = mean(West,2);
    end
    outeriter = outeriter + 1;
end
West = West * MAX; %re-scale to the original space

W_est = max(West,0);
H_est = FGMfcnls(X,W_est);
Y_res = W_est*H_est;


function Ax = afun(x, cond)
    % function for implicitly computing svds of X-v*e^T where X is a sparse
    % matrix. The function afun satisfies these required conditions:
    % Afun(x,'notransp') accepts a vector x and returns the product A*x.
    % Afun(x,'transp') accepts a vector x and returns the product A'*x.

    if strcmp(cond,'notransp')
        Ax = X * x - v * sum(x);
    else
        Ax = X' * x - v' * x;
    end
end


end


%%
function [Z_tilde,Theta,Delta,flag] = algorithm2_update_theta(Y,r,lambda0,Z_tilde,options)

if nargin <= 4
    options = [];
end
if ~isfield(options,'maxiter')
    options.maxiter = 5; 
end
if ~isfield(options,'timelimit')
    options.timelimit = 60; 
end
% initialization
n = size(Y,2);
p = r;
m = r-1;
Theta =zeros(r-1,r);
% initialize Z (and theta)
[Z_tilde,flag] = initial_theta(Y,r,lambda0,Z_tilde,options);
% checking whether initial optimization problem was successful
if ~flag
    Theta = randn(r-1,r);
    Delta = randn(r,n);
    return 
end
Z_0 = Z_tilde;
iter = 1;
Z_pre = randn(p,r);
Delta = [];
H = sparse(p+m+n+m,p+m+n+m);
H(p+m+1:end-m,p+m+1:end-m)=speye(n);
while iter < options.maxiter ... 
        && norm(Z_tilde - Z_pre,'fro')/norm(Z_pre,'fro')>1e-3
    Z_pre = Z_tilde;
    iter = iter + 1;
    for i = 1 : r %updating each column of Z sequentially
        lambda = lambda0 * det(Z_tilde)^2/det(Z_0)^2; %update lambda
        inv_Z = pinv(Z_tilde); % compute Z^-1
        zz = Z_tilde;
        zz(:,i)=[]; % taking the i-th (current) column out (for second equality constraint on alpha)
        f = inv_Z(i,:); %Z^-1(i,:), refer to the paper for the details
        f = [f zeros(1,m) zeros(1,n) zeros(1,m)]; %optimizaton parameters: [Z, Theta, Delta, alpha]
        A = [zeros(n,p) Y' -speye(n) zeros(n,m)]; %inequality constraint: Y' * Theta <= 1+ Delta --> Y' * Theta - Delta <= 1
        b = ones(n,1);
        Aeq1 = [eye(p) -eye(p,r-1) zeros(p,n) zeros(p,m)]; %equality constraint 1: Z = [Theta; e']
        beq1 = [zeros(m,1);1];
        Aeq2 = [zeros(m,p) eye(m) zeros(m,n), zz(1:m,:)]; %equality constraint 2: 0 to be in the interior of estimated Thetas
        beq2 = zeros(m,1);
        opts = optimoptions('quadprog','Display','none'); 
        [Z_tilde_i,~,e,~] = quadprog(lambda*H,-f',A,b,[Aeq1;Aeq2],[beq1;beq2],[-inf*ones(p+m+n,1); 0.01*ones(m,1)],[],[],opts);
        if e<0 && e~=-4 %check whether optimization was successful
            Z_tilde = randn(p,r);
            flag = 0;
            return
        else %update i-th column of Z, Theta, Delta and i-th element of alpha
            Z_tilde(:,i) = Z_tilde_i(1:p);
            theta = Z_tilde_i(p+1:p+m);
            Delta(:,i) = Z_tilde_i(p+m+1:end-m);
            Theta(:,i) = theta;
            alpha(:,i) = Z_tilde_i(end-m+1:end);
        end
    end
end
flag = 1;
end
%%
function [Z_tilde,flag] = initial_theta(Y,r,lambda0,Z_tilde,options)

%cputime0 = cputime; 
if nargin <= 4
    options = [];
end
if ~isfield(options,'maxiter')
    options.maxiter = 100; 
end
% if ~isfield(options,'timelimit')
%     options.timelimit = 60; 
% end
%initialization
n = size(Y,2);
p = r;
m = r-1;
Theta =zeros(r-1,r);
iter = 1;
Z_pre = randn(p,r);
Delta = [];
% initializing H matrix for the QP optimization problem
H = sparse(p+m+n+m,p+m+n+m);
H(p+m+1:end-m,p+m+1:end-m)=speye(n);
lambda = lambda0;
%main loop
while iter < options.maxiter ... 
        && norm(Z_tilde - Z_pre,'fro')/norm(Z_pre,'fro')>1e-2 %...
        %&& cputime-cputime0 <= options.timelimit 
    Z_pre = Z_tilde;
    iter = iter + 1;
    for i = 1 : r %updating each column of Z sequentially
        inv_Z = pinv(Z_tilde); %compute z^-1
        zz = Z_tilde;
        zz(:,i)=[]; % taking the i-th (current) column out
        f = inv_Z(i,:); %Z^-1(i,:), refer to the paper for the details
        f = [f zeros(1,m) zeros(1,n) zeros(1,m)]; %optimizaton parameters: [Z, Theta, Delta, alpha]
        A = [zeros(n,p) Y' -speye(n) zeros(n,m)]; %inequality constraint: Y' * Theta <= 1+ Delta --> Y' * Theta - Delta <= 1
        b = ones(n,1); 
        Aeq1 = [eye(p) -eye(p,r-1) zeros(p,n) zeros(p,m)]; %equality constraint 1: Z = [Theta; e']
        beq1 = [zeros(m,1);1];
        Aeq2 = [zeros(m,p) eye(m) zeros(m,n), zz(1:m,:)]; %equality constraint 2: 0 to be in the interior of estimated Thetas
        beq2 = zeros(m,1);
        opts = optimoptions('quadprog','Display','none'); %Turn off display information of the optimizer
        [Z_tilde_i,~,e,~] = quadprog(lambda*H,-f',A,b,[Aeq1;Aeq2],[beq1;beq2],[-inf*ones(p+m+n,1); 0.01*ones(m,1)],[],[],opts);
        if e<0 && e~=-4 % check if optimization has failed
            Z_tilde = zeros(p,r);
            flag = 0;
            return
        else %update i-th column of Z, Theta, Delta and i-th element of alpha
            Z_tilde(:,i) = Z_tilde_i(1:p);
            theta = Z_tilde_i(p+1:p+m);
            Delta(:,i) = Z_tilde_i(p+m+1:end-m);
            Theta(:,i) = theta;
            alpha(:,i) = Z_tilde_i(end-m+1:end);
        end
        
    end
end
flag = 1;
end
%%
function [c] = is_correct(Theta_c)
[m,num_c] = size(Theta_c);
f = [zeros(num_c,1)]; %no objective function, this is only a feasibility checking, parameter: c \in R^num_c
Aeq1 = [Theta_c]; % equality constraint 1: Theta_c * c = 0
beq1 = zeros(m,1);
Aeq2 = [ones(1,num_c)]; % equality constraint 2: c'*e = 1
beq2 =1;
Aeq = [Aeq1;Aeq2];
beq = [beq1;beq2];
lb = zeros(num_c,1); % c>=0 : convexity constraint
ub = [ones(num_c,1)]; % c<=1 : convexity constraint
[sol,fval] = linprog(-f,[],[],Aeq,beq,lb,ub); %LP with constraints
try
    c = sol(1:num_c); %if successful return coefficients
catch
    c = nan; %if not successful return nan
end
end

%%
function [V,e,UtU,UtM] = FGMfcnls(M,U,V,maxiter) 

% Fast gradient method to solve least squares on the unit simplex.  
% See Nesterov, Introductory Lectures on Convex Optimization: A Basic 
% Course, Kluwer Academic Publisher, 2004. 
% 
% This code solves: 
% 
%             min_{V(:,j) in Delta, forall j}  ||M-UV||_F^2, 
%
% where Delta = { x | sum x_i <= 1, x_i >= 0 for all i }.
%  
% See also Appendix A in N. Gillis, Successive Nonnegative Projection 
% Algorithm for Robust Nonnegative Blind Source Separation, arXiv, 2013. 
% 
%
% [V,e] = FGMfcnls(M,U,V,maxiter) 
%
% ****** Input ******
% M      : m-by-n data matrix
% U      : m-by-r basis matrix
% V      : initialization for the fast gradient method 
%          (optional, use [] if none)
% maxiter: maximum numbre of iterations (default = 500). 
%
% ****** Output ******
% V      : V(:,j) = argmin_{x in Delta}  ||M-Ux||_F^2 forall j. 
% e      : e(i) = error at the ith iteration

[m,n] = size(M); 
[m,r] = size(U); 
% Initialization of V 
if nargin <= 2 || isempty(V)
    V = zeros(r,n); 
    for i = 1 : n
        % Distance between ith column of M and columns of U
        disti = sum( (U - repmat(M(:,i),1,r)).^2 ); 
        [a,b] = min(disti); 
        V(b,i) = 1; 
    end
end
if nargin <= 3
    maxiter = 500; 
end
nM = norm(M,'fro')^2;
% Hessian and Lipschitz constant 
UtU = U'*U; 
L = norm(UtU,2); 
% Linear term 
UtM = U'*M; 
nM = norm(M,'fro')^2; 
alpha0 = 0.05; % Parameter, can be tuned. 
alpha(1) = alpha0;
V = SimplexProj( V ); % Project initialization onto the simplex
Y = V; % second sequence
i = 1; 
% Stop if ||V^{k}-V^{k+1}||_F <= delta * ||V^{0}-V^{1}||_F
delta = 1e-6;
eps0 = 0; eps = 1;  
while i <= maxiter && eps >= delta*eps0
    % Previous iterate
    Vp = V; 
    % FGM Coefficients  
    alpha(i+1) = ( sqrt(alpha(i)^4 + 4*alpha(i)^2 ) - alpha(i)^2) / (2); 
    beta(i) = alpha(i)*(1-alpha(i))/(alpha(i)^2+alpha(i+1)); 
    % Projected gradient step from Y
    V = SimplexProj( Y - (UtU*Y-UtM) / L );
    % `Optimal' linear combination of iterates
    Y = V + beta(i)*(V-Vp); 
    % Error 
    e(i) = nM - 2*sum(sum(V.*(UtM))) + sum(sum((UtU).*(V*V'))); 
    % Restart: fast gradient methods do not guarantee the objective
    % function to decrease, a good heursitic seems to restart whenever it
    % increases although the global convergence rate is lost! This could
    % be commented out. 
    if i >= 2 && e(i) > e(i-1)
        Y = V; 
    end 
    if i == 1
        eps0 = norm(V-Vp,'fro'); 
    end
    eps = norm(V-Vp,'fro'); 
    i = i + 1; 
end  
end
%%
function x = SimplexProj(y)

% Given y,  computes its projection x* onto the simplex 
% 
%       Delta = { x | x >= 0 and sum(x) <= 1 }, 
% 
% that is, x* = argmin_x ||x-y||_2  such that  x in Delta. 
% 
%  
% See Appendix A.1 in N. Gillis, Successive Nonnegative Projection 
% Algorithm for Robust Nonnegative Blind Source Separation, arXiv, 2013. 
% 
%
% x = SimplexProj(y)
%
% ****** Input ******
% y    : input vector.
%
% ****** Output ******
% x    : projection of y onto Delta.

x = max(y,0); 
K = find(sum(x) > 1); 
x(:,K) = blockSimplexProj(y(:,K));

end

%%
function x = blockSimplexProj(y)

% Same as function SimplexProj except that sum(max(Y,0)) > 1. 
[r,m] = size(y); 
ys = sort(-y);  ys = -ys;
indi2 = 1:m; lambda = zeros(1,m); 
S(1,indi2) = 0; 
for i = 2 : r
    if i == 2
        S(i,:) = (ys(1:i-1,:)-repmat(ys(i,:),i-1,1)); 
    else
        S(i,:) = sum(ys(1:i-1,:)-repmat(ys(i,:),i-1,1)); 
    end
    indi1 = find(S(i,indi2) >= 1); 
    indi1 = indi2(indi1);
    indi2 = find(S(i,:) < 1);
    if ~isempty(indi1)
        if i == 1
            lambda(indi1) = -ys(1,indi1)+1;
        else
            lambda(indi1) = (1-S(i-1,indi1))/(i-1) - ys(i-1,indi1);
        end
    end
    if i == r
        lambda(indi2) = (1-S(r,indi2))/r - ys(r,indi2);
    end
end
x = max( y + repmat(lambda,r,1), 0); 
end