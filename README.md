# Frisch Problem - Minimize the Rank of the Covariance Matrix Based on Noisy Data
For an example of the Frisch Problem, based on random "observations," run main.m.
To solve the Frisch Problem on an m x m, rank r matrix with T observations and random noise between 0 and 1 added (as Sigma_tilde) and using Sigma^-1 as the weight, run

	[rank_target, trace_min, noise, D] = TraceMinObsWeight( m, T, r, 1)

To solve the Frisch Problem on an m x m, rank r matrix with T observations and random noise between 0 and 1 added (as Sigma_tilde) and using an interative reweighting process, run
	
	[rank_target, trace_min, noise, D] = TraceMinObsWeightIter( m, T, r, 1)

Background on this is provided in the Analysis folder.
