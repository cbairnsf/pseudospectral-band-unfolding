%	terry_method_m_is_three
%This function implements the analytic version of the n --> \infty 1d SSH model
%Right now the relative scaling factor kappa is always 1. 
%
%Inputs:
%	v, w -- hopping amplitudes of the SSH model
%	kmin, kmax, Emin, Emax, npts -- specification of search space
%	eig_flag -- either the string "evals" or "evects"
%	folded_flag -- either the string "folded" or "unfolded"
%
%Outputs:
%	data -- a 3d array of sqrt_sigma_min values for heatmap plotting
%		if "evects" is selected then this is a 4d array such that (:,:,:,1) is sqrt_sigma_min values
%		and (:,:,:,2:end) are the associated eigenvectors
%
%Example usage:
%	data = terry_method_m_is_three(0.7, 1.4, -pi, pi, -4, 4, 101, "evals", "folded")
function data = terry_method_m_is_three(v, w, kmin, kmax, Emin, Emax, npts, eig_flag, folded_flag)
	%Test input and throw errors if it's out of bounds
	validate_input(v, w, kmin, kmax, Emin, Emax, npts, eig_flag, folded_flag);

%Doing these as inline functions in this way makes the code much slower for, presumably, computer reasons
%	H = @(k) [[0, v + exp(-i*k) * w]; [v + exp(i*k) * w, 0]]; 	%Hamiltonian (intra-cell part)
%	T = @(k) [[0, 1]; [exp(i*k), 0]];				%Translation op (intra-cell part)
%	Q = @(A, k) quadratic_gap({H(k), T(k)}, A, 1, "dense", "evals");%Quadratic gap as anonymous function
		%A is a length 2 vector, the (E, k) probe-site
		%k is the "internal" k
%We invert v and w to minimize the plot at k = 0, which is helpful to Alex for materials science reasons
	v = -v; 
	w = -w; 
	kappa = 1;

	%Loop over search space
	k_range = linspace(kmin, kmax, npts);	%Probe sites
	k_length = length(k_range);
	E_range = linspace(Emin, Emax, npts);
	E_length = length(E_range);
	%Majority of CPU time is spent in this parallel loop
	p = gcp('nocreate');	%Check for existing parpool
	if isempty(p)		%I should make this a separate m file since I use this code in multiple scripts
		[~, numprocs] = system('nproc');	%Two lines avoids shell output
		workers = str2num(numprocs) - 2;	%Parallel workers
		p = parpool("Processes", workers);
	end
	opts = parforOptions(p,'RangePartitionMethod','fixed','SubrangeSize',2);
	pctRunOnAll warning off;	%Disable warnings for parallel workers
	loop_tracking = tic;
	prog_bar_exists = maybe_create_progress_bar(k_length); 
	parfor (kndex = 1:k_length, opts)	
		%H = [[0, v + exp(-i*k_range(kndex)) * w]; [v + exp(i*k_range(kndex)) * w, 0]]; 
		%T = [[0, 1]; [exp(i*k_range(kndex)), 0]]; 
		H = [[0, v, w * exp(-i*k_range(kndex))]; [v, 0, v]; [w * exp(i*k_range(kndex)), v, 0]]; 
		T = [[0, 1, 0]; [0, 0, 1]; [exp(i*k_range(kndex)), 0, 0]]; 
		switch folded_flag
			case "folded"
				T = T * T * T; 	%Replace T with T^3 
			case "unfolded"
				%Do nothing
		end
		for jndex = 1:E_length
			M1 = (H - E_range(jndex) * eye(3))' * (H - E_range(jndex) * eye(3)) ;
			for index = 1:k_length
				M2 = (T - exp(i * k_range(index)) * eye(3))' * (T - exp(i * k_range(index)) * eye(3));
				Q = M1 + kappa * M2;
				switch eig_flag
					case "evals"
						data(jndex, index, kndex) = realsqrt(min(abs(eig(Q))));	%Not sparse because this is 2x2
					case "evects" %BROKEN RIGHT NOW
						[evects, evals] = eig(Q, "vector"); 
						[sigma_min, location] = min(evals);%Smallest eigenvalue and location
						best_eval = realsqrt(sigma_min);		%Smallest eigenvalue	
						best_evect = evects(:, location);	%Record best eigenvector
				end
				%temp = Q([E_range(jndex), k_range(index)], k_range(kndex));
				%data(jndex, index, kndex) = temp(1);	%Just keep eigenvalue
			end
		end
		if prog_bar_exists
			update_progress_bar; 
		end
	end
%	switch eig_flag
%		case "evals"
%			%do nothing
%		case "evects"
%			data(:,:,:,1) = eval_data;
%			data(:,:,:,2:3) = evect_data;
%	end
	str = sprintf("Total elapsed time: %f sec.", toc(loop_tracking));
	disp(str)
end

%Validate the input. Make sure values of the right type, in the allowable range, etc. 
function validate_input(v, w, kmin, kmax, Emin, Emax, step_size, eig_flag, folded_flag)
	err_code = "minimize_ssh_terry:validate_input";
	if ~(isscalar(v) && isscalar(w) && isscalar(kmin) && isscalar(kmax) &&...
		isscalar(Emin) && isscalar(Emax) && isscalar(step_size))
		err_msg = "All inputs must be scalars."
		error(err_code, err_msg)
	elseif ~(isnumeric(v) && isnumeric(w) && isnumeric(kmin) && isnumeric(kmax) &&...
		isnumeric(Emin) && isnumeric(Emax) && isnumeric(step_size))
		err_msg = "All inputs must be numeric."
		error(err_code, err_msg)
	elseif ~(isreal(v) && isreal(w) && isreal(kmin) && isreal(kmax) &&...
		isreal(Emin) && isreal(Emax) && isreal(step_size))
		err_msg = "All inputs must be real."
		error(err_code, err_msg)
	elseif ~(kmin <= kmax && Emin <= Emax && step_size > 0)
		err_msg = "Lower bounds must be less than upper bounds and step_size must be positive."
		error(err_code, err_msg)
	elseif ~(isstring(eig_flag))
		err_msg = "The last input must be either the string 'evals' or the string 'evects'."; 
	elseif ~(strcmp(eig_flag, "evals") || strcmp(eig_flag, "evects"))
		err_msg = "The last input must be either the string 'evals' or the string 'evects'."; 
	elseif ~(isstring(folded_flag))
		err_msg = "The last input must be either the string 'evals' or the string 'evects'."; 
	elseif ~(strcmp(folded_flag, "folded") || strcmp(folded_flag, "unfolded"))
		err_msg = "The last input must be either the string 'folded' or the string 'unfolded'."; 
	end
end

