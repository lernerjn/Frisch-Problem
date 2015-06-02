function [ rank_target, trace_min , noise , D] = TraceMinObsWeight( m, T ,r , rand_noise )
%Author: Jeremy Lerner, Stony Brook University
% TraceMinObsWeight estimates the minimum trace of a matrix, using
% the inverse of Sigma as the weight for the convex relaxation of
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
%       r: the desired rank of the matrix Sigma_hat
%       m: Sigma_hat will be an m x m matrix
%       range: the range for the random numbers in Sigma_hat
%       multiplier: a multiplier for the noise, Sigma_tilde
%
% Outputs:
%       rank_target: this should be r, the desired rank for Sigma_hat


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

[n,~] = size(Sigma);
% The first and only weight is (Sigma)^-1
W = inv(Sigma);

% Use CVX to optimize the trace, with variables Sigma_h (n x n matrix)
% and D (n x 1 diagonal matrix)
cvx_begin
    variables D(n) Sigma_h(n,n)
    minimize trace(W*(Sigma_h))
    subject to
        %This is equivalent to D>=0, as D is represented as a vector
        diag(D) == semidefinite(n);
        Sigma_h == semidefinite(n);
        Sigma == Sigma_h + diag(D);
cvx_end

trace_min = trace(W*Sigma_h);

end

