%	The one-dimensional SSH model
%
%This function computes the quadratic gap on a grid of points for the 1d generalized SSH model. 
%The region considered is [-pi, pi] x [-4, 4] in (k, E) ``probe site space''. 
%We use periodic boundary conditions.
%
%Inputs:
%	n -- the number of unit cells, each with m sites (positive integer)
%	m -- the number of sites per cell (positive integer)
%	v -- the intra-cell hopping amplitude (real number)
%	w -- the inter-cell hopping amplitude (real number)
%	kappa -- a parameter used to define the quadratic composite operator (real number)
%	npts -- the number of data points in each direction (positive integer s.t. the total number of points is npts^2)
%	pbc flag -- an optional input. Type "no pbc" to eliminate the periodic boundary conditions. 
%	T^m flag -- an optional input. Type "folded" to replace T with T^m. [H, T^m] = 0, but you get folded bands. 
%
%Output: 
%	H -- the Hamiltonian
%	T -- the translation operator 
%	data_out -- an npts x npts array of quadratic gaps
%
%Example usage:
%	[H, T, data_out] = ssh_model_improved(12, 4, 0.7, 1.4, 1, 101, "no pbc", "folded");
%
function [H, T, data_out] = ssh_model_improved(n, m, v, w, kappa, npts, varargin)
	%Test the inputs to make sure they are of correct type etc. 
	if ~(validate_input(n, "positive integer") && validate_input(m, "positive integer") && ...
		validate_input(npts, "positive integer") && validate_input(v, "real") && ...
		validate_input(w, "real") && validate_input(kappa, "real") && n > 0 && m > 1)	
		error("One or more input variables is not of an allowed type/value")
	end

	[pbc_flag, folded_flag] = test_input_strings(varargin);

	kmin = -pi; kmax = pi; Emin = -4; Emax = 4;	%Hardcoded bounds
	v = -v;		%This is useful to Alex for materials science reasons
	w = -w;		%For positive v, w, this centers the global minimum at k = 0

	%Now we create the Hamiltonian and translation operator
	row_indices = 1:(n * m);					%Row indices for the entries of H and T
	column_indices = [2:(n * m), 1];				%Column indices
	T = sparse(row_indices, column_indices, ones(1, n * m));	%The T matrix
	if ~pbc_flag		%Eliminate PBC term if requested
		T(n*m, 1) = 0;
	end
	if folded_flag		%This happens if any input was provided beyond the required six numbers
		T = T^m;	%[H, T^m] = 0--we should see a folded band in this case
	end
	sub_and_sup_diag = [repmat([repmat(v, 1, m - 1), w], 1, n - 1), repmat(v, 1, m - 1), w];	
	H = sparse(row_indices, column_indices, sub_and_sup_diag);	%The H matrix (before symmetrization)
	if ~pbc_flag		%Eliminate PBC term if requested
		H(n * m, 1) = 0;
		H(1, n * m) = 0;
	end
	H = H + H';							%Symmetrize it

	%The probe site coordinates to loop over and a container for data
	kvals = exp(i * linspace(kmin, kmax, npts));	%Map points in R to points in S^1 using exp
	Evals = linspace(Emin, Emax, npts);
	data_out = zeros(npts, npts);			%Need to allocate in advance so it's a "sliced variable"

	[~, opts] = maybe_create_parpool;		%Either use already existing parpool or create a new one	
	prog_bar = maybe_create_progress_bar(npts);	%Create a progress bar only in the desktop environment
	global_tracking_id = tic;			%Track total runtime in loop

	%====================the main data collection loop=====================
	parfor (index = 1:npts, opts)
		M2 = T - kvals(index) * speye(n * m);		%The matrix expressing the eigenvalue problem for T
		S2 = M2' * M2;					%Force it to be hermitian
		for jndex = 1:npts
			M1 = H - Evals(jndex) * speye(n * m);	%Same for H
			S1 = M1' * M1;				%Hermitian again
			Q = S1 + kappa^2 * S2;			%The quadratic composite operator
			eigval = eigs(Q, 1, 'smallestabs');	%Quadratic gap
			data_out(jndex, index) = realsqrt(real(eigval)) ;	%NB -- the order of indices is reversed
		end
		if prog_bar
			update_progress_bar();			%Increment progress bar if it exists
		end
	end
	%======================================================================

	str = sprintf("Total elapsed time: %f", toc(global_tracking_id));	%Show total runtime
	disp(str)
end


function [pbc_flag, folded_flag] = test_input_strings(string_array)
	pbc_flag = 1;						%Defaults
	folded_flag = 0;
	num_strings = size(string_array, 2); 			%Number of varargin inputs
	for index = 1:num_strings				%Test each one for being a string
		if ~isstring(string_array{index})
			error("Too many non-string inputs.")
		else
			switch string_array{index}		%For strings, test them vs options
				case "no pbc"
					pbc_flag = 0;
				case "folded"
					folded_flag = 1;
				otherwise
					error("Unrecognized option.")
			end
		end
	end
end
