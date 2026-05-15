%=======spatially_resolved_breathing_heatmap============================================================================
%
%This function does a loop over purported eigenvalues of H, T1, T2, and T3 from the breathing honeycomb model. 
%(Or any model, but that's the intent.) Create the operators with the c6v code from Alex. 
%
%Inputs:
%	H, T1, T2, T3, X, Y	--	operators from c6v code (we don't need chiral mirror things)
%	kappa_T, kappa_XY	--	kappa parameters for Q (kappa_H is 1 always) (must be real numbers)
%	pos			--	a 1x2 vector or real numbers indicating the position for spatial resolution 
%	npts 			--	number of data points in each direction (a total of npts^3) (positive integer)
%	sr_flag			--	spatially resolved flag; 0 means don't do spatial resolution, 1 means do do it
%
%Output: 
%	data_out 		-- 	an npts x npts x npts array of quadratic gaps
%
%Example usage:
%	data_out = spatially_resolved_breathing_heatmap(H, T1, T2, T3, X, Y, kappa_X, kappa_XY, pos, npts, sr_flag)
%	data_out = spatially_resolved_breathing_heatmap(H, T1, T2, T3, X, Y, 2.0, 0.5, [12.1, 13.5], 128, 1)
%
function data_out = spatially_resolved_breathing_heatmap(H, T1, T2, T3, X, Y, kappa_T, kappa_XY, pos, npts, sr_flag)
	check_inputs(H, T1, T2, T3, X, Y, kappa_T, kappa_XY, pos, npts, sr_flag)	%Make sure inputs are all okay

	%Hard-coded global variables (a = 1 needed to match c6v code)
	kxmin = -pi; kymin = -pi; kxmax =  pi; kymax =  pi; Emin = -4; Emax = 4; a = 1; 

	%Create variables
	E_list = linspace(Emin, Emax, npts);	%Energies to use
	kx_list = linspace(kxmin, kxmax, npts);	%Angles to use; they become points later 
	ky_list = linspace(kymin, kymax, npts);	%Angles to use; they become points later 
	a1 = a * [3/2,  sqrt(3)/2];		%Fundamental vectors of the regular hexagon
	a2 = a * [3/2, -sqrt(3)/2]; 
	a3 = a1 + a2; 				%Needs to be + because c6v code does T3 = T2 * T1; 
	data_out = zeros(npts, npts, npts);	%Preallocate so it's a "sliced" variable
	num_sites = size(H, 1);			%Number of sites
	MX = X - pos(1) * speye(num_sites); 	%Eigenvalue problem for X (doesn't change in loop so done here)
	MY = Y - pos(2) * speye(num_sites); 	%And for Y (same comment)
	MXYs = sr_flag * kappa_XY' * kappa_XY * (MX' * MX + MY' * MY);	%Done once, not thousands of times in loop

	%Set up loop, progress indicator, and run-time tracking
	[p, opts] = maybe_create_parpool; 	%Use already existing parpool or create one and get parpool opts
	prog_bar = maybe_create_progress_bar(npts);			%Create progress bar only in desktop environment
	global_tracking_id = tic;		%Using output variable allows nesting tic's

	%======================Main Data Collection Loop================================================================
	parfor (index = 1:npts, opts)
		for jndex = 1:npts
			kxky_vec = [kx_list(index), ky_list(jndex)]	%Location in pc
			T1_eval = exp(i * kxky_vec * a1'); 	%Exp map [-pi, pi] --> S^1 turns angles into points
			M1 = T1 - T1_eval * speye(num_sites); 	%Construct eigenvalue problem from purported eigenvalue
			T2_eval = exp(i * kxky_vec * a2'); 
			M2 = T2 - T2_eval * speye(num_sites);	%Eigenvalue problem for T2
			T3_eval = exp(i * kxky_vec * a3');	
			M3 = T3 - T3_eval * speye(num_sites);	%Eigenvalue problem for T3
			for kndex = 1:npts
				M4 = H - E_list(kndex) * speye(num_sites);%Eigenvalue problem for H
				%Make the quadratic composite operator
				Q = M4' * M4 + kappa_T' * kappa_T * (M1' * M1 + M2' * M2 + M3' * M3) + MXYs;	
				eigval = eigs(Q, 1, "smallestabs");	%Grab the quadratic gap
				data_out(jndex, index, kndex) = realsqrt(real(eigval));	%NB -- order of indices!
			end
		end
		if prog_bar
			update_progress_bar();			%Increment progress bar if it exists
		end
	end
	%===============================================================================================================

	str = sprintf("Total elapsed time: %f", toc(global_tracking_id));	%Display runtime
	disp(str)
end
%=======================================================================================================================

%=======check_inputs====================================================================================================
%This function checks the inputs to make sure they are of correct type etc. 
function check_inputs(H, T1, T2, T3, X, Y, kappa_T, kappa_XY, pos, npts, sr_flag)
	validateattributes(H, {'numeric'}, {'2d', 'square', 'nonempty'})	%A non-empty square matrix of numbers 
	size_H = size(H);							%Grab the size of H for comparing
	validateattributes(T1, {'numeric'}, {'2d', 'square', 'size', size_H})	%A square matrix the same size as H
	validateattributes(T2, {'numeric'}, {'2d', 'square', 'size', size_H})	%A square matrix the same size as H
	validateattributes(T3, {'numeric'}, {'2d', 'square', 'size', size_H})	%A square matrix the same size as H
	validateattributes(X, {'numeric'}, {'2d', 'square', 'size', size_H})	%A square matrix the same asize as H & T
	validateattributes(Y, {'numeric'}, {'2d', 'square', 'size', size_H})	%A square matrix the same asize as H & T
	validateattributes(kappa_T, {'numeric'}, {'scalar', 'real'})		%Real number
	validateattributes(kappa_XY, {'numeric'}, {'scalar', 'real'})		%Real number
	validateattributes(pos, {'numeric'}, {'vector', 'real', 'size', [1, 2]})%A vector
	validateattributes(npts, {'numeric'}, {'scalar', 'integer', 'positive'})%A positive integer
	validateattributes(sr_flag, {'numeric'}, {'scalar', 'integer'})		%A one entry character array 
	assert(sr_flag == 0 || sr_flag == 1)					%Flag is either 0 or 1
end	
%=======================================================================================================================

