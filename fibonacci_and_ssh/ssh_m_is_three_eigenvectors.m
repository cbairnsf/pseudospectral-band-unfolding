%	ssh_m_is_three_eigenvectors
%This function implements m = 3 trivalent SSH model and extracts eigenvector data along the analytically unfolded curve.
%Right now the following values are hard-coded:
%	kappa = 1
%	v = 0.8, w = 1.2
%	m = 3
%	k ranges from -pi to pi
%	We are using PBCs
%
%Inputs:
%	N 	- the number of unit *cells* (so 3N sites) (positive integer)
%	npts	- the number of points to consider anlong the anaytically unfolded curve (positive multiple of 3)
%
%Outputs:
%	coords	- a npts x 2 length array of k values respresenting the probe space coords used as (k, E) pairs
%	H	- the Hamiltonian 
%	T	- the translation operator
%	evects 	- an N x npts array where data(:, i) is the eigenvector at the ith point considered
%	gaps	- a 1 x npts array where gaps(i) is the quadratic gap at the ith point considered
%
%Example usage:
%	[coords, H, T, evects, gaps] = ssh_m_is_three_eigenvectors(24, 512)
function [coords, H, T, evects, gaps] = ssh_m_is_three_eigenvectors(N, npts)
	check_input(N, npts);							%Make sure inputs are allowable
	
	kappa = 1;								%Hardcoded for now
	v = 0.8;
	w = 1.2;
	m = 3; 
	k_vals = linspace(-pi, pi, npts);					%k values to consider
	k_evals = exp(i * k_vals);						%the corresponding eigenvalues
	%[E_vals, ~] = ssh_model_exact_solution(v, w, m, k_vals); 		%E values to consider, exact
	E_vals = abs(v + w) * cos(k_vals); 					%E values to consider, approximate
	coords = [k_evals.', E_vals.']; 					%Coordinates of probe sites
	num_sites = m * N;							%Total number of sites
	evects = zeros(num_sites, npts);					%Initialize eigenvector output
	gaps = zeros(1, npts);							%Initialize gaps output

	row_indices = 1:num_sites; 						%Data for sparse T matrix
	column_indices = [2:num_sites, 1]; 
	T = sparse(row_indices, column_indices, ones(1, num_sites));		%The T matrix
	sub_and_sup_diag = [repmat([v, v, w], 1, N - 1), v, v, w];		%Data for sparse H matrix
	H = sparse(row_indices, column_indices, sub_and_sup_diag);		%The H matrix (before symmetrization)
	H = H + H'; 								%Symmetrize

	[parp, opts] = maybe_create_parpool; 					%Use existing parpool or create one 
	prog_bar = maybe_create_progress_bar(npts);				%Create a progress bar in the DE only
	global_tracking_id = tic;						%Track runtime as a diagnostic

	%=======main data collection loop===============================================================================
	parfor (index = 1:npts, opts)						%The loop
		k = k_evals(index);						%Purported eigenvalues
		E = E_vals(index);
		M1 = H - E * speye(num_sites); 					%Eigenvalue problem for H
		M1s = M1' * M1; 						%Make it Hermitian
		M2 = T - k * speye(num_sites);					%Eigenvalue problem for T
		M2s = kappa^2 * M2' * M2; 					%Make it Hermitian and scale by kappa
		Q = M1s + M2s; 							%Quadratic composite operator
		[evec_Q, eval_Q] = alt_eigs(Q, 1); 				%Get eigenvalue and eigenvector
		evects(:, index) = evec_Q; 					%Save eigenvector
		gaps(index) = realsqrt(real(eval_Q));				%Save gap
		if prog_bar 
			update_progress_bar();					%Increment progress bar if it exists
		end
	end
	%===============================================================================================================

	str = sprintf("Total elapsed time: %f", toc(global_tracking_id));	%Display runtime
	disp(str)
end
%=======================================================================================================================

%=======check_inputs====================================================================================================
%Makes sure that input values are of allowable types and fall within acceptable ranges and extract number of sites
function check_input(N, npts)	
	validateattributes(N, {'numeric'}, {'scalar', 'integer', 'positive'})	%Positive integer
	validateattributes(npts, {'numeric'}, {'scalar', 'integer', 'positive'})%Positive integer
	assert((npts / 3) == floor(npts / 3)); 					%Divisible by three
end
%=======================================================================================================================

