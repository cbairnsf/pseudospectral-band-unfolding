%	alt_eigs
%
%This function calls eigs in a ``try-catch'' mode then, if the output was a NaN, calls it again with larger values for
%some of the parameters to try to fix that. It specifically fixes eigval = eigs(__, 1, 'sa') since that is the mode of
%operation used in all of the QCO code. This code is intended to be run hundreds of thousands of times, so there is no
%input checking at all. It's all about that speed. 
%
%	Input:	
%		in		--	an nxn sparse matrix
%		vec_flag 	--	0 or 1, to either discard or compute evects (not checked--don't input anything else!)
%	
%	Output:
%		evect		--	0 or else the eigenvector corresponding to the below eigenvalue
%		evalue 		--	the smallest eigenvalue of that matrix
%
%	Example usage:
%		[evect, evalue] = alt_eigs(in, 1);
%		[~, eval] = alt_eigs(in, 0);
%
function [evect, evalue] = alt_eigs(in, vec_flag)
	%Option structure starting values
	opts = struct; 
	opts.tol = 1e-6; 			%This is the default, maybe consider lowering it if runtime is an issue
	opts.maxit = 300;			%Iterations before failure
	opts.disp = 0;				%Write nothing to output
	%opts.issym = 0;			%Hermitian is not usually symmetric, and T isn't either anyway (might count as symmetric anyway from CHATGPT)
	num_dims = size(in, 1);			%How large is the state space?
	opts.p = min(num_dims - 1, 20);		%Default when we are only looking for one eigenvector

	success = 0;				%Have we found the eigenvector yet?
	evect = 0; 				%Initialize


	%The ``try-catch'' loop (this is not true try/catch, which does exist in MATLAB)
	if vec_flag
		while ~success 
			[evect, evalue] = eigs(in, 1, 'smallestabs', opts);	%Run eigs using the parameters
			if isnan(evalue)			%NaN indicates failure--double things as a result, then try again.
				%disp("It happened!")
				opts.maxit = 2 * opts.maxit;
				opts.p = min(2 * opts.p, num_dims - 1);	%Cannot exceed size of state space
			else
				success = 1;
			end
		end
	else
		while ~success 
			[~, evalue] = eigs(in, 1, 'smallestabs', opts);	%Run eigs using the parameters
			if isnan(evalue)			%NaN indicates failure--double things as a result, then try again.
				%disp("It happened!")
				opts.maxit = 2 * opts.maxit;
				opts.p = min(2 * opts.p, num_dims - 1);	%Cannot exceed size of state space
			else
				success = 1;
			end
		end
	end
end
