%=======smoothed_HTXQ_looping===========================================================================================
%This function accepts three operators, H, T, and X, as well as vector lists of 
%purported eigenvalues for these operators, and three ``kappa'' parameters, then
%constructs the quadratic composite operator
%	Q = k_H (H - E I)^2 + k_T (T - e^(ika) I)^2 + k_X f((X - x I))^2
%where ``^2'' means "hermitian square" when the operator in question is not 
%hermitian. Here f is a smooth bump function implemented in smooth_bumo_function.m. 
%Then all combinations of elements from the lists are substituted for
%E, k, and x and a 3d array of the square root of the smallest eigenvalue of the
%quadratic composite operator is produced. Note that the middle operator is 
%assumed to be unitary, with the list of purported eigenvalues being provided as
%complex arguments. Also note that a is assumed to be 1. 
%There are commented out lines, prefaced with '%%%%%', that switch the function 
%from outputing only qudratic gap values, or also outputting the six other values
%necessary to create the "peanut plots". 
%	Required external m-files:
%		maybe_create_parpool
%		maybe_create_progress_bar 	(in DE only)
%		update_progress_bar 		(in DE only)
%	Inputs:
%		H	--	a square numeric array
%		T	--	a square numeric array the same size as H
%		X	--	a square numeric array the same size as H and T (runs faster if NOT sparse?)
%		k_H	--	a real number
%		k_T	--	a real number
%		k_X	--	a real number
%		H_pev	--	a numeric vector
%		T_pev	--	a numeric vector (can be a different length)
%		X_pev	--	a numeric vector (can be a different length)
%	Ouputs:
%		data_out --	a length(H_pev) x length(T_pev) x length(X_pev)
%				numeric array of non-negative real numbers OR
%		data_out --	that x 7 where the other six entries are the 
%				expected values and variances of the three 
%				operators
%	Example Usage:
%		data_out = smoothed_HTXQ_looping(H, T, X, 1, 1, 0.5, Evals, kvals, xvals); 
function data_out = HTXQ_looping(H, T, X, k_H, k_T, k_X, H_pev, T_pev, X_pev)
	check_inputs(H, T, X, k_H, k_T, k_X, H_pev, T_pev, X_pev)		%Make sure inputs are of OK types etc.
	num_sites = size(H, 1);							%Total number of sites
	lH = length(H_pev); lT = length(T_pev); lX = length(X_pev);		%Lengths of each list of purp. eigvals
	total_num_iterations = lH * lT * lX;					%The total number of quadratic gaps 
	%X_smoothing_matrix = diag(sparse(smooth_bump_function(30, 34, 38, 42, diag(X), "inverse"))); %To smooth X later
	%X_smoothing_matrix = diag(sparse(smooth_bump_function(26, 31, 39, 46, diag(X), "inverse"))); %To smooth X later
	%X_smoothing_matrix = diag(sparse(smooth_bump_function(16, 26, 46, 56, diag(X), "inverse"))); %To smooth X later
	%X_smoothing_matrix = diag(sparse(smooth_bump_function(36, 72, 144, 180, diag(X), "inverse"))); %To smooth X later
	%X_smoothing_matrix = diag(sparse(smooth_bump_function(100, 125, 250, 275, diag(X), "inverse"))); %To smooth X later
	X_smoothing_diagonal = smooth_bump_function(900, 987, 987*2, 987*2 + 87, diag(X), "inverse"); %To smooth X later
	X_diagonal = full(diag(X)); 
	%X_smoothing_matrix = diag(ones(1, size(H, 1))); 
	%%%%%data_out = zeros(lH, lT, lX, 7);					%Initialize data_out to all zeros
	data_out = zeros(lH, lT, lX);						%Initialize data_out to all zeros
	[parp, opts] = maybe_create_parpool; 					%Use existing parpool or create one 
	prog_bar = maybe_create_progress_bar(num_sites);			%Create a progress bar in the DE only
	global_tracking_id = tic;						%Track runtime as a diagnostic
	%=======main loop===============================================================================================
	parfor (index = 1:lH, opts)						%Outer loop (E direction)
		if ~prog_bar
			loop_tracking_id = tic; 
		end
		E = H_pev(index);						%Purported eigenvalue of H
		M1 = (H - E * speye(num_sites));				%Matrix representing eigenvalue problem
		M1s = k_H^2 * M1' * M1;						%Hermitian square of the matrix w/ kappa
		for jndex = 1:lT						%Middle loop (k direction)
			k = exp(i * T_pev(jndex));				%Purported eigenvalue of T
			M2 = (T - k * speye(num_sites));			%Matrix representing eigenvalue problem
			M2s = k_T^2 * M2' * M2;					%Hermitian square with kappa parameter
			for kndex = 1:lX					%Inner loop (x direction)
				x = X_pev(kndex);				%Purported eigenvalue of X
				%M3 = (X - x * speye(num_sites));			%Matrix representing eigenvalue problem
				%M3 = X_smoothing_matrix * M3; 
				%M3s = k_X^2 * M3' * M3;				%Hermitian square with kappa parameter
				%M3s = k_X^2 * X_smoothing_matrix * M3' * M3;
				M3_diagonal = X_smoothing_diagonal .* (X_diagonal - x).^2; 
				M3s = k_X^2 * spdiags(M3_diagonal, 0, num_sites, num_sites); 
				Q = M1s + M2s + M3s;				%Create the quadratic composite operator
	%%%%%			[evect_Q, eval_Q] = eigs(Q, 1, 'sa');		%Compute the quadratic gap
				%eval_Q = eigs(Q, 1, 'smallestabs');		%Compute the quadratic gap
				[~, eval_Q] = alt_eigs(Q, 0); 			%Compute the quadratic gap, avoiding NaN 
				%expH = expect(H, evect_Q);			%Compute expectation values of all three
				%expT = expect(T, evect_Q);
				%expX = expect(X, evect_Q);
				%varH = expect(H' * H, evect_Q) - expH' * expH; 	%Compute variances of all three 
				%varT = expect(T' * T, evect_Q) - expT' * expT; 
				%varX = expect(X' * X, evect_Q) - expX' * expX; 
				%data_out(index, jndex, kndex, :) = [eval_Q;...	%Store the data
			%				expH - E;...
			%				expT - k;...
			%				expX - x;...
			%				varH;
			%				varT;
			%				varX];
				data_out(index, jndex, kndex) =...
						realsqrt(real(eval_Q));		%Add data to data_out using linear index
			end
		end
		if prog_bar 
			update_progress_bar();			%Increment progress bar if it exists
		else
			str = sprintf("Loop elapsed time: %f", toc(loop_tracking_id));	%Display runtime
			disp(str)
		end
	end
	%===============================================================================================================
	str = sprintf("Total elapsed time: %f", toc(global_tracking_id));	%Display runtime
	disp(str)
end
%=======================================================================================================================

%=======check_inputs====================================================================================================
%Makes sure that input values are of allowable types and fall within acceptable ranges
function check_inputs(H, T, X, k_H, k_T, k_X, H_pev, T_pev, X_pev)
	validateattributes(H, {'numeric'}, {'2d', 'square', 'nonempty'})	%A non-empty square matrix of numbers 
	size_H = size(H);							%Grab the size of H for comparing
	validateattributes(T, {'numeric'}, {'2d', 'square', 'size', size_H})	%A square matrix the same size as H
	validateattributes(X, {'numeric'}, {'2d', 'square', 'size', size_H})	%A square matrix the same asize as H & T
	validateattributes(k_H, {'numeric'}, {'scalar', 'real'})		%Real number
	validateattributes(k_T, {'numeric'}, {'scalar', 'real'})		%Real number
	validateattributes(k_X, {'numeric'}, {'scalar', 'real'})		%Real number
	validateattributes(H_pev, {'numeric'}, {'vector'})			%Vector of numbers
	validateattributes(T_pev, {'numeric'}, {'vector'})			%Vector of numbers
	validateattributes(X_pev, {'numeric'}, {'vector'})			%Vector of numbers
end
%=======================================================================================================================

%=======expect==========================================================================================================
%Find the expectation value of an operator in a state
function exptval = expect(M, v)
	exptval = (M' * v)' * v;						%This is the definition of expect. value
end
%=======================================================================================================================

