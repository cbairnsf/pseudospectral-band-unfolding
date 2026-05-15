function data_out = breathing_heatmap(H, T1, T2, T3, kappa, npts)
	%Create needed variables
	kxmin = -pi; kymin = -pi; kxmax =  pi; kymax =  pi; Emin = -4; Emax = 4;
	a = 3;
	E_list = linspace(Emin, Emax, npts);	%Energies to use
	kx_list = linspace(kxmin, kxmax, npts);	%Angles to use; they become points next 
	ky_list = linspace(kymin, kymax, npts);	%Angles to use; they become points next 
	a1 = [3/2,  sqrt(3)/2];			%Fundamental vectors of the regular hexagon
	a2 = [3/2, -sqrt(3)/2]; 
	a3 = a1 - a2; 
	data_out = zeros(npts, npts, npts);	%Preallocate so it's a "sliced" variable
	num_sites = size(H, 1);			%Number of sites--this slightly depends on angle

	[p, opts] = maybe_create_parpool; 	%Use already existing parpool or create one and get parpool opts
	prog_bar = maybe_create_progress_bar(npts);	%Create a progress bar only in the desktop environment
	global_tracking_id = tic;		%Using output variable allows nesting tic's

	%======================Main Data Collection Loop=======================
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
				M4 = H - E_list(kndex) * speye(num_sites);	%Eigenvalue problem for H
				Q = M4' * M4 + kappa' * kappa * (M1' * M1 + M2' * M2 + M3' * M3);%Quadratic composite operator
				eigval = real(eigs(Q, 1, "smallestabs"));	%Grab the quadratic gap
				data_out(jndex, index, kndex) = realsqrt(output_data); 	%NB -- order of indices!
			end
		end
		if prog_bar
			update_progress_bar();			%Increment progress bar if it exists
		end
	end
	%======================================================================

	str = sprintf("Total elapsed time: %f", toc(global_tracking_id));	%Display runtime
	disp(str)
end
