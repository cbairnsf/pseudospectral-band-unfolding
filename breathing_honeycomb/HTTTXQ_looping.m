%=======HTTTXQ_looping==================================================================================================
%This function accepts six operators, H, T1, T2, T3, X, and Y, as well as lists of 
%purported eigenvalues for these operators, and three ``kappa'' parameters, then
%constructs the quadratic composite operator. (The kappas for all three T's are equal.)
%	Q = k_H (H - E I)^2 + \sum_i k_T_i (T_i - e^(i k_i a) I)^2 + k_X (X - x I)^2
%where ``^2'' means "hermitian square" when the operator in question is not 
%hermitian. Then all combinations of elements from the lists are substituted for
%E, k_i, and x and a 3d array of the square root of the smallest eigenvalue of the
%quadratic composite operator is produced. Note that the middle three operators are 
%assumed to be unitary, with the list of purported eigenvalues being provided as
%complex arguments. Also note that a is assumed to be 1. 
%	Required external m-files:
%		maybe_create_parpool
%		maybe_create_progress_bar 	(in DE only)
%		update_progress_bar 		(in DE only)
%		alt_eigs
%	Inputs:
%		H	--	a square numeric array
%		T1	--	a square numeric array the same size as H
%		T2	--	a square numeric array the same size as H
%		T3	--	a square numeric array the same size as H
%		X	--	a square numeric array the same size as H and T
%		Y	--	a square numeric array the same size as H and T
%		k_H	--	a real number
%		k_T	--	a real number
%		k_X	--	a real number
%		k_Y	--	a real number
%		H_pev	--	a numeric vector
%		T_pev	--	a numeric vector (can be a different length, is used for both T1 and T2)
%		X_pev	--	a numeric vector (can be a different length)
%		Y_pev	--	a numeric vector (can be a different length)
%	Ouputs:
%		data_out --	a length(T_pev) x length(T_pev) x length(H_pev) x length(X_pev) x length(Y_pev)
%				numeric array of non-negative real numbers 
%	Example Usage:
%		data_out = HTTTXQ_looping(H, T1, T2, T3, X, Y, 1, 1, 0.5, 0.5, Evals, kvals, kvals, xvals, yvals); 
function data_out = HTTTXQ_looping(H, T1, T2, T3, X, Y, k_H, k_T, k_X, k_Y, H_pev, T_pev, X_pev, Y_pev)
	check_inputs(H, T1, T2, T3, X, Y, k_H, k_T, k_X, k_Y,...		%Make sure inputs are of OK types etc.
						H_pev, T_pev, X_pev, Y_pev)
	num_sites = size(H, 1);							%Total number of sites
	%This is the balanced triangular approach, not the original approach!
	a_1 = -[3/2, sqrt(3)/2];						%Basis vector for a real-space unit cell
	a_2 = [3/2, -sqrt(3)/2];						%Basis vector for a real-space unit cell
	a_3 = -a_1 - a_2;							%This vector is also helpful
	lH = length(H_pev); lT = length(T_pev);					%Lengths of each list of purp. eigvals
	lX = length(X_pev); lY = length(Y_pev);		
	total_num_iterations = lH * lT * lT * lX * lY;				%The total number of quadratic gaps 
	data_out = zeros(lT, lT, lH, lX, lY);					%Initialize data_out to all zeros
	[parp, opts] = maybe_create_parpool; 					%Use existing parpool or create one 
	prog_bar = maybe_create_progress_bar(total_num_iterations);		%Create a progress bar in the DE only
	global_tracking_id = tic;						%Track runtime as a diagnostic
	%=======main loop===============================================================================================
	%parfor (index = 1:lH, opts)						%Outer loop (E direction)
	for index = 1:lH							%Switch parfor where best 
		E = H_pev(index);						%Purported eigenvalue of H
		M1 = (H - E * speye(num_sites));				%Matrix representing eigenvalue problem
		M1s = k_H^2 * M1' * M1;						%Hermitian square of the matrix w/ kappa
		%for jndex = 1:1T
		parfor (jndex = 1:lT, opts) 					%Middle loop (k direction)
			for kndex = 1:lT
				k = [T_pev(jndex), T_pev(kndex)];		%A point in 2d k-space
				kev = [exp(i * k * a_1'),...			%Related to purported eigenvalue of T's
					exp(i * k * a_2'), exp(i * k * a_3')];			
				M2 = (T1 - kev(1) * speye(num_sites));		%Matrix representing eigenvalue problem
				M2s = k_T^2 * M2' * M2;				%Hermitian square with kappa parameter
				M3 = (T2 - kev(2) * speye(num_sites));		%Matrix representing eigenvalue problem
				M3s = k_T^2 * M3' * M3;				%Hermitian square with kappa parameter
				M4 = (T3 - kev(3) * speye(num_sites));		%Matrix representing eigenvalue problem
				M4s = k_T^2 * M4' * M4;				%Hermitian square with kappa parameter
				for lndex = 1:lX				%Inner loop (x direction)
					x = X_pev(lndex);			%Purported eigenvalue of X
					M5 = (X - x * speye(num_sites));	%Matrix representing eigenvalue problem
					M5s = k_X^2 * M5' * M5;			%Hermitian square with kappa parameter
					for mndex = 1:lY
						y = Y_pev(mndex);		%Purported eigenvalue of X
						M6 = (Y - y * speye(num_sites));%Matrix representing eigenvalue problem
						M6s = k_Y^2 * M6' * M6;		%Hermitian square with kappa parameter
						Q = M1s + M2s + M3s +...	%Create the quadratic composite operator
							M4s + M5s + M6s;
						%[~, eval_Q] = alt_eigs(Q, 0);	%Compute the quadratic gap
						eval_Q = eigs(Q, 1, 'smallestabs');
						data_out(kndex, jndex,...	%Note the order!
						index, lndex, mndex) =...	%Remove +0i parts
						realsqrt(real(eval_Q));		%Take sqrt
						if prog_bar 
							update_progress_bar();	%Increment progress bar if it exists
						end
					end
				end
			end
		end
	end
	%===============================================================================================================
	str = sprintf("Total elapsed time: %f", toc(global_tracking_id));	%Display runtime
	disp(str)
end
%=======================================================================================================================

%=======check_inputs====================================================================================================
%Makes sure that input values are of allowable types and fall within acceptable ranges
function check_inputs(H, T1, T2, T3, X, Y, k_H, k_T, k_X, k_Y, H_pev, T_pev, X_pev, Y_pev)
	validateattributes(H, {'numeric'}, {'2d', 'square', 'nonempty'})	%A non-empty square matrix of numbers 
	size_H = size(H);							%Grab the size of H for comparing
	validateattributes(T1, {'numeric'}, {'2d', 'square', 'size', size_H})	%A square matrix the same size as H
	validateattributes(T2, {'numeric'}, {'2d', 'square', 'size', size_H})	%A square matrix the same size as H
	validateattributes(T3, {'numeric'}, {'2d', 'square', 'size', size_H})	%A square matrix the same size as H
	validateattributes(X, {'numeric'}, {'2d', 'square', 'size', size_H})	%A square matrix the same asize as H & T
	validateattributes(Y, {'numeric'}, {'2d', 'square', 'size', size_H})	%A square matrix the same asize as H & T
	validateattributes(k_H, {'numeric'}, {'scalar', 'real'})		%Real number
	validateattributes(k_T, {'numeric'}, {'scalar', 'real'})		%Real number
	validateattributes(k_X, {'numeric'}, {'scalar', 'real'})		%Real number
	validateattributes(k_Y, {'numeric'}, {'scalar', 'real'})		%Real number
	validateattributes(H_pev, {'numeric'}, {'vector'})			%Vector of numbers
	validateattributes(T_pev, {'numeric'}, {'vector'})			%Vector of numbers
	validateattributes(X_pev, {'numeric'}, {'vector'})			%Vector of numbers
	validateattributes(Y_pev, {'numeric'}, {'vector'})			%Vector of numbers
end
%=======================================================================================================================
