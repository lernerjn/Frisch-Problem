%Author: Jeremy Lerner, Stony Brook University
% Note: this file takes about 4 minutes to run. The speed can be
% drastically increased by using a parfor on line 39 of
% TraceMinObsTrials.m, for which this file takes about 2 minutes

% The entries in this row vector are the sizes of matrices that will be
% used in the rank minimization programs.
n_in = [10, 100];

% The number of observations used to create the noise free covariance 
% matrix, as in, using 1.1*n observations.
T = 1.1;

%Call the file that averages the results for 
[ average_inverse, average_iter ] = TraceMinObsTrials( T , n_in );

j = 1;
s = cell(1,2*numel(n_in));
for i = 1:2*numel(n_in)
    if mod(i,2) == 0
        s(i) = {sprintf('Actual_Rank_for_%0.5g',n_in(j))}; 
        j = j+1;
    else
        s(i) = {sprintf('Trace_for_%0.5g',n_in(j))}; 

    end
end

% Display the average output for using the inverse as the weight
Average_Using_Simple_Inverse_Table = array2table(average_inverse','VariableNames', s)

% Display the average output for using the iterative weight method
Average_Using_Iterative_Weight_Table = array2table(average_iter','VariableNames', s)
