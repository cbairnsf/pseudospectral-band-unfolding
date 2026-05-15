%	fibonacci_two_values
%
%This function computes the quadratic gap on a grid of points for the 1d generalized SSH model for a quasi-crystal 
%based on the Fibonacci sequence. 
%The region considered is [-pi, pi] x [-4, 4] in (k, E) ``probe site space''. 
%We use periodic boundary conditions.
%
%Inputs:
%	n -- the number of unit cells, each with F_m sites (positive integer, should probably be 1)
%	m -- number of iterations on the starting word, i.e. length(word) = F_m
%	repl -- the "replacement word" used to define the recursion (format 'vwwvwwwvw'; see fibonacci_word help)
%	v -- the A to B hopping amplitude (real number)
%	w -- the B to A hopping amplitude (real number)
%	kappa -- a parameter used to define the quadratic composite operator (real number)
%	npts -- the number of data points in each direction (positive integer s.t. the total number of points is npts^2)
%	pbc_flag -- an optional input. Type "no pbc" to enforce open boundary conditions. 
%	folded_flag -- an optional input. Type "folded" to replace T with T^m. [H, T^m] = 0, but you get folded bands. 
%	ops_flag -- an optional input. Type "just operators" to just get H and T without computing a heatmap at all. 
%	provided_flag -- an optional input. Type "operators provided" to read in H and T rather than recomputing. 
%		if you use this option, you must input H and T as the final two inputs
%
%Output: 
%	H -- the Hamiltonian
%	T -- the translation operator 
%	data_out -- an npts x npts array of quadratic gaps
%
%Example usage:
%	[H, T, data_out] = fibonacci_two_values(1, 12, 'LLS', 0.8, 1.2, 0.1, 101);
%
%Terry's Instructions: 
%	Start with the word LLS 
%	Take n = 1 (no periodicity)
%
function [H, T, data_out] = ssh_model_fibonacci_terry(n, m, repl, v, w, kappa, npts, varargin)
	%Test the inputs to make sure they are of correct type etc. 
	check_inputs(n, m, repl, v, w, kappa, npts)

	%Set flags based on input strings (if any)
	[pbc_flag, folded_flag, ops_flag, provided_flag] = test_input_strings(varargin);

	kmin = -pi; kmax = pi; Emin = -2; Emax = 2;			%Hardcoded bounds
	v = -v;		%This is useful to Alex for materials science reasons
	w = -w;		%For positive v, L, S, this centers the global minimum at k = 0

	%Now we create the Hamiltonian and translation operator
	word = fibonacci_word(m, 'v', 'w', repl); 			%Fibonacci word for n cells (external m file)
	if ~provided_flag
		row_indices = 1:(n * length(word));				%Row indices for the entries of H and 
		column_indices = [2:(n * length(word)), 1];			%Column indices
		T = sparse(row_indices, column_indices, ones(1, n * length(word)));	%The T matrix
		if ~pbc_flag							%Eliminate PBC term if requested
			T(n * length(word), 1) = 0;
		end
		if folded_flag	%This happens if any input was provided beyond the required six numbers
			%T = T^(length(word));%[H, T^m] = 0--we should see a folded band in this case
			T = T^2;
		end
		list_of_commas = repmat(',', 1, length(word));			
		shuffle = [word; list_of_commas];   %Need commas and v's in the word 
		decorated_word = shuffle(1:(end - 1));				%Of the form L,S,L,S,S,L,S,...
		array_word = ['[', decorated_word, ']']; 			%You can tell this will be bad...
		one_unit_cell_pattern = eval(array_word);			%Remember, eval is evil. Use with care. 
		sub_and_sup_diag = repmat(one_unit_cell_pattern, 1, n); 	%Repeat n times
		H = sparse(row_indices, column_indices, sub_and_sup_diag);	%The H matrix (before symmetrization)
		H = H + H';	%Symmetrize it
		if ~pbc_flag	%Eliminate PBC term if requested
			H(n * length(word), 1) = 0;
			H(1, n * length(word)) = 0;
		end
	else	%In this case the H and T operators need to be the last two inputs in varargin 
		assert(length(varargin) >= 2)
		[H, T] = extract_ops(varargin(end-1:end));
	end

	%The probe site coordinates to loop over and a container for data
	kvals = exp(i * linspace(kmin, kmax, npts));	%Map points in R to points in S^1 using exp
	Evals = linspace(Emin, Emax, npts);
	data_out = zeros(npts, npts);			%Need to allocate in advance so it's a "sliced variable"

	if ~ops_flag
		[~, opts] = maybe_create_parpool;		%Either use already existing parpool or create a new one	
		prog_bar = maybe_create_progress_bar(npts);	%Create a progress bar only in the desktop environment
		global_tracking_id = tic;			%Track total runtime in loop
	end

	%====================the main data collection loop=====================
	if ~ops_flag
		parfor (index = 1:npts, opts)
			M2 = T - kvals(index) * speye(n * length(word));%The eigenvalue problem for T
			S2 = M2' * M2;					%Force it to be hermitian
			for jndex = 1:npts
				M1 = H - Evals(jndex) * speye(n * length(word));	%Same for H
				S1 = M1' * M1;				%Hermitian again
				Q = S1 + kappa^2 * S2;			%The quadratic composite operator
				eigval = eigs(Q, 1, 'smallestabs');	%Quadratic gap
				%[eigvect, eigval] = alt_eigs(in, 1);
				%[~, eigval] = alt_eigs(Q, 0);
				data_out(jndex, index) = realsqrt(real(eigval));	%NB -- the order of indices is reversed
			end
			if prog_bar
				update_progress_bar();			%Increment progress bar if it exists
			end
		end
		str = sprintf("Total elapsed time: %f", toc(global_tracking_id));	%Show total runtime
		disp(str)
	else
		%Do nothing--output is already zero
	end
	%======================================================================

end

%This function makes sure the inputs are valid
function check_inputs(n, m, repl, v, w, kappa, npts)
	validateattributes(n, {'numeric'}, {'scalar', 'positive', 'integer'})	%A positive integer
	validateattributes(m, {'numeric'}, {'scalar', 'positive', 'integer'})	%A positive integer
	validateattributes(repl, {'char'}, {'scalartext', 'nonempty'})		%A one entry character array 
	validateattributes(v, {'numeric'}, {'scalar', 'real'})			%A real number
	validateattributes(w, {'numeric'}, {'scalar', 'real'})			%A real number
	validateattributes(kappa, {'numeric'}, {'scalar', 'real'})		%A real number
	validateattributes(npts, {'numeric'}, {'scalar', 'positive', 'integer'})%A positive integer
	for index = 1:length(repl)
		if (repl(index) ~= 'v') && (repl(index) ~= 'w') 
			error("Extraneous character in replacement word.")
		end
	end
end

%This function sets the flags based on however many strings were input
function [pbc_flag, folded_flag, ops_flag, provided_flag] = test_input_strings(string_array)
	pbc_flag = 1;						%Defaults
	folded_flag = 0;
	ops_flag = 0; 
	provided_flag = 0; 
	num_strings = size(string_array, 2); 			%Number of varargin inputs
	for index = 1:num_strings				%Test each one for being a string
		test_string = isstring(string_array{index}); 
		if ~test_string && ~provided_flag
			disp(string_array{index})
			error("Too many non-string inputs.")
		elseif test_string
			switch string_array{index}		%For strings, test them vs options
				case "no pbc"
					pbc_flag = 0;
				case "folded"
					folded_flag = 1;
				case "just operators"
					ops_flag = 1; 
				case "operators provided"
					provided_flag = 1;
				otherwise
					error("Unrecognized option.")
			end
		end
	end
end

%This function extracts H and T from varargin and makes sure they are of the correct type etc. 
function [H, T] = extract_ops(input_array)
	H = input_array{1}; 
	T = input_array{2}; 
	validateattributes(H, {'numeric'}, {'2d', 'square'}); 
	validateattributes(T, {'numeric'}, {'2d', 'square'}); 
	assert(size(H, 1) == size(T, 1))
end

