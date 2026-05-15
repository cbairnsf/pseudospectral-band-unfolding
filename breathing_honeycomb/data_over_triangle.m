%	data_over_triangle
%
%This function computes quadratic gap values for given operators over a specified triangular prism in probe space. In 
%theory you can use any four matrices (of the same size) but the intended use case is a Hamiltonian and three 
%translation operators corresponding to a graphene-adjacent model (breathing, TBG, etc). 
%
%	Inputs:
%		H	--	the first operator (a square numeric array that should be hermitian)
%		Ti	--	the various translation operators (square numeric arrays of the same size, i = 1, 2, 3)
%		kappa_H	--	the kappa value used to form Q associtated with H (a non-negative real number)
%		kappa_T	--	the kappa value used to form Q associtated with Ti (i = 1, 2, 3) (also real >= 0)
%		Emin	--	the minimum purported eigenvalue of H to consider (a real number)
%		Emax	--	the maximum (a real number which is not less than than Emin)
%		G	--	a point in kx-ky space (a 2x1 numeric array of real numbers)
%		M	--	''
%		K	--	''
%		npts_E	--	the number of points above each point of the journey in the E direction (positive int)
%		npts_k	--	the number of points on the triangular path in k-space
%							(positive int, slightly modified if necessary to fit triangle)
%
%	Outputs:
%		data_out	--	an npts_E x 3 * npts_k array of quadratic gap values
%	
%	Example usage:
%		data_out = data_over_triangle(H, T1, T2, T3, kappa_H, kappa_T, Emin, Emax, G, M, K, npts_E, npts_k)
%
%The points for a Gamma-->M-->K-->Gamma triangle in the reciprocal honeycomb lattice (for a = 1) are
%	G = [0; 0];								%The Gamma point is the origin
%	M = [2 * pi / 3;  0];							%Location of a hexagon edge midpoint
%	K = [2 * pi / 3; 2 * pi / (3 * sqrt(3))];				%Location of a Dirac point
%The actual path is a 30-60-90 right triangle. 
%
function data_out = data_over_triangle(H, T1, T2, T3, kappa_H, kappa_T, Emin, Emax, G, M, K, npts_E, npts_k)
	num_sites = check_inputs(H, T1, T2, T3, kappa_H, kappa_T,...
					Emin, Emax, G, M, K, npts_E, npts_k);	%Validate input and get size of system

	%Create the k_vals for looping and the a vectors
	[k1, k2, k3] = scaled_triangle(npts_k);					%Get best possible side lengths
	npts_k = k1 + k2 + k3;							%Change user input as needed
	%Get sidelengths (in pixels)
	first_leg = [linspace(G(1), M(1), k1); linspace(G(2), M(2), k1)];	%First side of triangle
	second_leg = [linspace(M(1), K(1), k2); linspace(M(2), K(2), k2)];	%Second side of triangle
	third_leg = [linspace(K(1), G(1), k3); linspace(K(2), G(2), k3)];	%Third side of triangle
	k_vals = [first_leg, second_leg, third_leg]; 				%The whole triangle in kxky-space
	E_vals = linspace(Emin, Emax, npts_E);					%The purported energy eigenvalues
	jvals = 1:npts_E;							%Dumb extra variable because ``sliced''
	a_1 = -[3/2, sqrt(3)/2];						%Basis vector for a real-space unit cell
	a_2 = [3/2, -sqrt(3)/2];						%Basis vector for a real-space unit cell
	a_3 = -a_1 - a_2;							%This vector is also helpful

	%Initialize output and tracking variables before data collection
	data_out = zeros(npts_E, npts_k);					%Initialize data_out to all zeros
	[parp, opts] = maybe_create_parpool; 					%Use existing parpool or create one 
	prog_bar = maybe_create_progress_bar(npts_k);				%Create a progress bar in the DE only
	global_tracking_id = tic;						%Track runtime as a diagnostic

	%=======main data collection loop===============================================================================
	parfor (index = 1:npts_k, opts)						%Outer loop (k direction)
	%for index = 1:npts_k
		k = k_vals(:, index);						%Location in kxky-space
		kev = [exp(i * a_1 * k), exp(i * a_2 * k),...
		       				exp(i * a_3 * k)];		%Related to purported eigenvalue of T's
		M1 = (T1 - kev(1) * speye(num_sites));				%Matrix representing eigenvalue problem
		M1s = kappa_T^2 * M1' * M1;					%Hermitian square with kappa parameter
		M2 = (T2 - kev(2) * speye(num_sites));				%Matrix representing eigenvalue problem
		M2s = kappa_T^2 * M2' * M2;					%Hermitian square with kappa parameter
		M3 = (T3 - kev(3) * speye(num_sites));				%Matrix representing eigenvalue problem
		M3s = kappa_T^2 * M3' * M3;					%Hermitian square with kappa parameter
		MT = M1s + M2s + M3s;					
		for jndex = jvals						%Inner loop (E direction)
			E = E_vals(jndex);					%Purported eigenvalue of H
			M4 = (H - E * speye(num_sites));			%Matrix representing eigenvalue problem
			MH = kappa_H^2 * M4' * M4;				%Hermitian square of the matrix w/ kappa
			Q = MH + MT;						%Create the quadratic composite operator
			%eval_Q = eigs(Q, 1, 'sa');				%Compute the quadratic gap
			[~, eval_Q] = alt_eigs(Q, 0); 
			data_out(jndex, index) = real(eval_Q);			%Add data to data_out
		end
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
function num_sites = check_inputs(H, T1, T2, T3, kappa_H, kappa_T, Emin, Emax, G, M, K, npts_E, npts_k)	
	validateattributes(H, {'numeric'}, {'2d', 'square', 'nonempty'})	%A non-empty square matrix of numbers 
	assert(all(all(H == H')))						%Check hermiticity
	size_H = size(H);							%Grab the size of H for comparing
	num_sites = size_H(1);							%Output total number of sites to main
	validateattributes(T1, {'numeric'}, {'2d', 'square', 'size', size_H})	%A square matrix the same size as H
	validateattributes(T2, {'numeric'}, {'2d', 'square', 'size', size_H})	%A square matrix the same size as H
	validateattributes(T3, {'numeric'}, {'2d', 'square', 'size', size_H})	%A square matrix the same size as H
	if ~(issparse(H) && issparse(T1) && issparse(T2) && issparse(T3))	%Warn against using full matrices
		disp("Using non-sparse arrays may result in long runtime!")
	end
	validateattributes(kappa_H, {'numeric'},...
	       				{'scalar', 'real', 'nonnegative'})	%Real number >= 0
	validateattributes(kappa_T, {'numeric'},...
	       				{'scalar', 'real', 'nonnegative'})	%Real number >= 0
	validateattributes(Emin, {'numeric'}, {'scalar', 'real'})		%A real number
	validateattributes(Emax, {'numeric'}, {'scalar', 'real'})		%A real number
	assert(Emax >= Emin)							%Need Emin <= Emax
	validateattributes(G, {'numeric'}, {'vector', 'real', 'size', [2, 1]})	%Need G, M, and K to be 2x1 real vectors
	validateattributes(M, {'numeric'}, {'vector', 'real', 'size', [2, 1]})	%Need G, M, and K to be 2x1 real vectors
	validateattributes(K, {'numeric'}, {'vector', 'real', 'size', [2, 1]})	%Need G, M, and K to be 2x1 real vectors
	validateattributes(npts_E, {'numeric'},...
					{'scalar', 'integer', 'positive'})	%Positive integer
	validateattributes(npts_k, {'numeric'},...
	       				{'scalar', 'integer', 'positive'})	%Positive integer
end
%=======================================================================================================================

%=======scaled_triangle=================================================================================================
%Takes in an integer, the perimeter of a triangle, and finds the closest perimeter of a 30-60 right triangle to input.
%Outputs the side lengths of that triangle: k1 = long leg, k2 = short leg, k3 = hypotenuse. 
function [k1, k2, k3] = scaled_triangle(npts_k) 
	basic_perimeter = 3 + sqrt(3); 						%Perimeter of 1-2-sqrt(3) triangle
	scale_factor = round(npts_k / basic_perimeter, 0); 			%Nearest integer to a multiple of perim.
	k1 = round(sqrt(3) * scale_factor, 0); 					%Long leg
	k2 = scale_factor; 							%Short leg
	k3 = 2 * scale_factor;							%Hypotenuse
end
