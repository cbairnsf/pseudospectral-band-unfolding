%=======================================================================================================================
%	idos
%
%This script computes the density of states and integrated density of states for a supplied Hamiltonian.
%	Inputs:
%		H	--	the Hamiltonian, a Hermitian matrix
%		npts_E	--	number of energies to check between min(spec(H)) and max(spec(H)), a positive integer
%		sigma	--	parameter in gaussian
%
%	Outputs:
%		E_vals	--	the E values considered (for plotting)
%		dos	--	the desity of states
%		intdos	--	the integrated density of states
%		
%	Example Usage:
%		[E_vals, dos, intdos] = idos(H, 4001, 0.05)	%Values from alex calc2d script
%
function [E_vals, dos, intdos] = idos(H, npts_E, sigma)
	check_inputs(H, npts_E, sigma)						%Ensure that inputs are valid

	evals = real(eig(full(H)));						%Get all eigenvalues	
	E_min = min(evals) - 0.5; E_max = max(evals) + 0.5;			%Get min and max for energy axis
	E_vals = linspace(E_min, E_max, npts_E);				%All E values
	dos = zeros(size(E_vals));						%Initialize
	intdos = zeros(size(E_vals));						%Initialize

	for index = 1:length(evals)
        	dos = dos + exp(-(E_vals - evals(index)).^2 / (2 * sigma^2));	%Add Gaussian to graph
	end

	intdos(1) = dos(1);							%They start equal
	for jndex = 2:length(dos)
		intdos(jndex) = intdos(jndex - 1) + dos(jndex);			%Accumulation
	end
end
%=======================================================================================================================

%=======================================================================================================================
%	check_inputs
%This function makes sure the inputs are valid
%
function check_inputs(H, npts_E, sigma)
	validateattributes(H, {'numeric'}, {'2d', 'square'})			%A square matrix of numeric values
	Hdag = H';								%Hermitian conjugate
	assert(all(Hdag(:) == H(:))); 						%Insist of Hermtian to ensure real spec
	validateattributes(npts_E, {'numeric'}, {'scalar', 'integer',...	%A positive integer 
							'positive'})	
	validateattributes(sigma, {'numeric'}, {'scalar', 'real', 'positive'})	%A positive real number
end
%=======================================================================================================================

