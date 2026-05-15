%	peanut_plot_data_from_eigenvectors
%This function compute the five things plotted in a ``peanut plot'' when the gaps and corresponding eigenvectors are 
%already known. 
%
%Inputs:
%	H		- the Hamiltonian
%	T		- the translation operator
%	coords		- an npts x 2 vector of probe site coordinates, (k, E) pairs
%	evects		- a num_sites * npts array where evects(:, i) is the eigenvector at the ith probe site
%	gaps		- a 1 * npts array where gaps(i) is the quadratic gap at the ith probe site
%
%Outputs:
%	peanut_data	- a 5 x npts array where peanut_data(:, i) is all five relevant quantities at the ith probe site
%	
%The order of peanut_data(:, i) is:
%		gaps(i)
%		deviation of expectation of H in state evects(:, i) from purported eigenvalue E = coords(i, 2)
%		deviation of expectation of T in state evects(:, i) from purported eigenvalue k = coords(i, 1)
%		variance of H in state evects(:, i)
%		variance of T in state evects(:, i)
%
%The ``peanut corollary'' says that: 
%	peanut_data(1)^2 = peanut_data(2)^2 + kappa^2 * peanut_data(3)^2 + peanut_data(4) + kappa^2 * peanut_data(5) 
%	you need to remember what kappa is!
%
%Example urage:
%	peanut_data = peanut_plot_data_from_eigenvectors(H, T, coords, evects, gaps)
%
function peanut_data = peanut_plot_data_from_eigenvectors(H, T, coords, evects, gaps)
	npts = check_input(H, T, coords, evects, gaps);				%Make sure inputs are okay; get npts
	peanut_data = zeros(5, npts);						%Initialize

	[~, opts] = maybe_create_parpool;					%Make a parpool unless one exists
	prog_bar_exists = maybe_create_progress_bar(npts);			%Progress bar in the desktop env.
	global_tracking_id = tic;						%Track total runtime in loop

	%======================Main Data Collection Loop================================================================
	parfor (index = 1:npts, opts)
		evect_Q = evects(:, index); 					%Best state at probe site 
		k = coords(index, 1); 						%Probe site coords (note order!)
		E = coords(index, 2);		
		expH = expectation(H, evect_Q);					%Compute the needed things
		expT = expectation(T, evect_Q);
		varH = expectation(H' * H, evect_Q) - expH' * expH;
		varT = expectation(T' * T, evect_Q) - expT' * expT; 
		peanut_data(:, index) = [0;...					%Store most of the data, not gaps yet
					expH - E;...				%Deviations of exptected values from 
					expT - k;...				%probe site coordinates
					varH;...				%Variances
					varT];
		if prog_bar_exists
			update_progress_bar();					%Increment progress bar if it exists
		end
	end
	peanut_data(1, :) = gaps; 						%Include gaps in output
	%===============================================================================================================

	%Display information
	str = sprintf("Total elapsed time: %f", toc(global_tracking_id));		
	disp(str)							%Display runtime
end
%=======================================================================================================================

%=======check_inputs====================================================================================================
%Makes sure that input values are of allowable types and fall within acceptable ranges and extract number of probe sites
function length_gaps = check_input(H, T, coords, evects, gaps)	
	validateattributes(H, {'numeric'}, {'2d', 'square', 'nonempty'})	%A non-empty square matrix of numbers 
	size_H = size(H);							%Grab the size of H for comparing
	validateattributes(T, {'numeric'}, {'2d', 'square', 'size', size_H})	%A square matrix the same size as H
	validateattributes(gaps, {'numeric'}, {'vector', 'real'})		%Vector of numbers
	length_gaps = size(gaps, 2); 						%Grab the length of gaps for comparing
	validateattributes(evects, {'numeric'}, {'2d', 'size',...		%Numeric array of correct size
						[size_H(1), length_gaps]})
	validateattributes(coords, {'numeric'}, {'2d', 'size',...		%Numeric array of correct size
						[length_gaps, 2]})
end
%=======================================================================================================================

%======expectation======================================================================================================
%Finds the expectation value of an operator in a state
function exptval = expectation(M, v)
	exptval = v' * M * v;
end
%=======================================================================================================================
