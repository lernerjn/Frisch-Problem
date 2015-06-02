# Minimum Rank Matrix Completion
For an example of minimum rank matrix completion, based on a random n x n rank r matrix (where 3/4 of the entries are sampled), run

	[X, rankX, M, Q] = Min_Rank_Matrix_Completion( n, n, 3*n^2/4, 1e-6)

To run this code on an existing matrix, M, where the non-sampled entries are NaN, run
	
	Note: Q = [ 1 1 1; 2 2 2; 3 3 3] is the matrix we are trying to recover.
	M = [ nan nan 1; 2 2 nan ; 3 3 3];
	[X, rankX, M] = Min_Rank_Matrix_Completion( M, 1e-6);

Background on this is provided in the Analysis folder.
