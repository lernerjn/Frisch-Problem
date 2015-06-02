function [ rank_target, trace_min , noise , D] = TraceMinObsWeightIter( m, T ,r , rand_noise )
%Author: Jeremy Lerner, Stony Brook University
% TraceMinObsWeightIter estimates the minimum trace of a matrix, using an
% iterative weighting method for the convex relaxation of
% the rank minimization problem:
%
%   mr_+(Sigma) = min{ rank(Sigma_hat) | Sigma = Sigma_tilde + Sigma_hat, 
%                       Sigma_tilde and Sigma_hat are positive definite 
%                       and Sigma_tilde is diagonal}
%
%   convex relaxation, trace minimization:
%   min{ trace(Sigma - D) | Sigma, D and Sigma - D are positive definite,
%           and D is diagonal}
%
% Inputs:
%       
%       m: Sigma_hat will be an m x m matrix
%       T: the number of observations 
%       r: the desired rank of the matrix Sigma_hat
%       rand_multiplier: a multiplier for the noise, Sigma_tilde
%
% Outputs:
%       rank_target: the desired rank for Sigma_hat, the optimal point for
%           the convex optimization problem
%       trace_min: the actual estimate of the rank
%       noise: Sigma_tilde, as computed before the convex optimization, for
%           reference
%       D: the D used to find the optimal trace
%       


% Create the observations, using a normal distribution
x = zeros(m,T);
for i = 1:T
    x(:,i) = randn(m,1);
end

%Create r unique columns
for i = r+1:m
    x(i,:) = x(i-r,:);
end

%Create a rank r matrix, Sigma_hat
Sigma_hat = x*x';

%The optimal value for the convex relaxation is the rank of this matrix
%Sigma_hat
rank_target = rank(Sigma_hat);

% Sigma_tilde is the covariance matrix of the noise (and is diagonal)
noise = rand_noise*diag(rand(m,1));
Sigma = Sigma_hat + noise;
epsilon = 5e-3;

[n,~] = size(Sigma);
% Initial weight
W = inv(Sigma+epsilon*eye(n,n));
W_old = zeros(n,n);
D_old = ones(n,1);
D = zeros(n,1);
trace_min_old = 0;
trace_min = 1;

%Stop when either W, D or the estimate stops changing significnatly, and
% update W iteratively
while ( norm(W_old - W) > 1e-2 && norm(D_old - D) > 1e-2 && abs(trace_min - trace_min_old) > 1e-3 )
    W_old = W;
    D_old = D;
    trace_min_old = trace_min;
    
    % Use CVX to optimize the trace, with variables Sigma_h (n x n matrix)
    % and D (n x 1 diagonal matrix)
    cvx_begin
        variables D(n) Sigma_h(n,n)
        minimize trace(W*Sigma_h)
        subject to
            %This is equivalent to D>=0, as D is represented as a vector
            diag(D) == semidefinite(n);
            Sigma_h == semidefinite(n);
            Sigma == Sigma_h + diag(D);
    cvx_end
    trace_min =  trace(W*Sigma_h);
    W = inv(Sigma - diag(D) + epsilon*eye(n,n));
    
end


end

