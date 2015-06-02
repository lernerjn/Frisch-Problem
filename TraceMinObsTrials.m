function [ average_inverse, average_iter ] = TraceMinObsTrials( T , n_in )
%Author: Jeremy Lerner, Stony Brook University
% TraceMinObsTrials uses the functions TraceMinObsWeight and 
% TraceMinObsWeightIter to estimate the rank minimization problem:
%
%    mr_+(Sigma) = min{ rank(Sigma_hat) | Sigma = Sigma_tilde + Sigma_hat, 
%                       Sigma_tilde and Sigma_hat are positive definite 
%                       and Sigma_tilde is diagonal}
%
% Inputs:   
%           T: the percentage of n, where Sigma is an n x n matrix,
%           dictating how many 
%           
%           n_in: a vector of sizes to test the algorithm on matrices of
%           size n_in(i) x n_in(i) for i=1..numel(n_in). Note: the input
%           ranks increase as r = ceil(n*(j/10)), where j=1..10, that is,
%           they increase 10% each time, so for numbers not divisible by
%           10, the output is not as clear. Also, the algorithms work best
%           for matrices of size 10 x 10 or larger
%
%
% Outputs:
%           average_inverse: a matrix with the average of three trials for
%           the rank minimization problem of various sizes, using the
%           inverse as the weight
%
%           average_iter: a matrix with the average of three trials for
%           the rank minimization problem of various sizes, using an
%           iterative process to find an optimal weight

    i = -1;
    average_inverse = zeros(2*numel(n_in),10);
    average_iter = zeros(2*numel(n_in),10);
%     average_time = zeros(4,10);
    for n = n_in
        i = i + 2;
        
        % Optional: can use a parfor loop here instead to speed up the runtime
        parfor j=1:10
            % T_n is the number of observations, T_n = n*T
            T_n = round(n*T);
            
            % r is the rank of the noise free covariance matrix 
            r = ceil(n*(j/10));
            
            % Perform 3 trials, average the results
            for k = 1 :3
%                 tic
                % Call the file that runs the trace minimization using the
                % simple inverse as a weight
                [ ~, trace_min , ~ , ~] = TraceMinObsWeight( n, T_n ,r , 1 );
                average_inverse(i,j) = ((k-1)*average_inverse(i,j) + trace_min) / (k);
%                 average_time(i,j) = ((k-1)*average(i,j) + toc) / (k);
                
                % Call the file that runs the trace minimization using the
                % iterative method to find the best weight
                [ ~, trace_min , ~ , ~] = TraceMinObsWeightIter( n, T_n ,r , 1 );
                average_iter(i,j) = ((k-1)*average_iter(i,j) + trace_min) / (k);


            end
        end
    end
    
    %Fill in the rest of the table with the actual ranks
    i = 0;
    for n = n_in
        i = i + 2;
        for j=1:10
            r = ceil(n*(j/10));
            for k = 1 :3
                average_inverse(i,j) = r;
                average_iter(i,j) = r;
            end
        end
    end
    
end
