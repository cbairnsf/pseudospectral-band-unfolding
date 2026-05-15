%	ssh_model_exact_solution
%=======This function calculates the exact unfolded dispersion curve for the 1d SSH model===============================
%Using exact solutions from https://opg.optica.org/oe/fulltext.cfm?uri=oe-29-26-42827&id=465762 in the m = 3 case. 
%Citation for source for exact solutions in m = 3 case. 
%	@article{Zhang:21,
%		author = {Yiqi Zhang and Boquan Ren and Yongdong Li and Fangwei Ye},
%		journal = {Opt. Express},
%		keywords = {Femtosecond lasers; Numerical simulation; Phase; Quantum information; Reflection; Refractive index},
%		number = {26},
%		pages = {42827--42836},
%		publisher = {Optica Publishing Group},
%		title = {Topological states in the super-SSH model},
%		volume = {29},
%		month = {Dec},
%		year = {2021},
%		url = {https://opg.optica.org/oe/abstract.cfm?URI=oe-29-26-42827},
%		doi = {10.1364/OE.445301},
%	}
%
%	Inputs:
%		v, w	--	the hopping amplitudes
%		m 	--	the number of sites per unit cell (must be 2 or 3)
%		k 	--	the independent variable, a row vector of points in the PBZ (to be modded to [-pi, pi])
%				(NB when m = 3 length(k) must be divisible by 3)
%	
%	Outputs:
%		E	--	a vector of values--the unfolded dispersion curve
%		E_cell	--	a cell array of vectors of values--the folded dispersion curves
%	
%	Example Usage:
%		[E, E_cell] = ssh_model_exact_solution(0.7, 1.4, 3, linspace(-pi, pi, 128));
function [E, E_cell] = ssh_model_exact_solution(v, w, m, k)
	check_inputs(v, w, m, k) 			%Test inputs to ensure they are allowable

	switch m 					%Switch to cover various options for m
		case 2
			%We use different branches of the analytic dispersion curve on the right/left half of the circle
			%We use 2k here because we are mapping S^1 to itself via z --> z^2, i.e. k --> 2k
			k = k + pi;					%Centers minimum for materials science reasons
			k_wound_twice = imag(log(exp(i*k)));		%Map everything to the interval [-pi, pi] 
			E_pos = sqrt(v^2 + w^2 + 2*v*w * cos(2 * k));	%The two branches of the folded curve	
			E_cell = {E_pos, -E_pos}; 			%Combine into a cell array for output
			E = zeros(1, size(k, 2)); 			%Initialize output to correct length
			right_half = abs(k_wound_twice) < pi/2;		%Right half of circle (as a set of indices)
			E(right_half) = E_pos(right_half); 		%The positive branch on the right half
			E(~right_half) = -E_pos(~right_half);		%The negative branch on the left half
		case 3
			%Totally different technique--found analytic solution in a paper
			assert(length(k) / 3 == floor(length(k) / 3))	%Make sure that length k is divisible by 3
			k = k + pi;					%Centers minimum for materials science reasons
			k = imag(log(exp(i*k)));			%Map everything to the interval [-pi, pi] 
			p = (2 * v^2 + w^2) / 3;			%Auxilliary variables to simplify expressions
			q = v^2 * w * cos(k); 				%Here is the k dependence
			c2 = (q + sqrt(q.^2 - p^3)).^(1/3);
			c1 = p ./ c2; 					%c1 and c2 are coefficients needed below
			z1 = exp(i * pi/3);				%Complex numbers
			z2 = exp(-i * pi/3);				
			e1 = real(c1 + c2); 				%The three solutions
			e2 = real(-z1 * c1 - z2 * c2); 			%They should be real, but MATLAB sees 0i and 
			e3 = real(-z2 * c1 - z1 * c2); 			%treats it as complex. 
			E_cell = {e1, e2, e3}; 				%Combine them into a cell array for output
			E_unfolded = [e1, e2, e3]; 
			E_unfolded = E_unfolded(1:3:end); 
			%disp([size(E_unfolded); size(e1); size(e2); size(e3); size(E(left_third)); size(E(middle_third)); size(E(right_third)); size(E)])
			left_third = (k <= -pi/3); 
			middle_third = ((k > -pi/3) & k < (pi/3));
			right_third = (k >= pi/3);
			E = zeros(1, length(k)); 
			E(left_third) = E_unfolded(left_third); 
			E(middle_third) = E_unfolded(middle_third);
			E(right_third) = E_unfolded(right_third); 
		otherwise
			str = sprintf("The code has not been expanded to cover the case m = %d yet, sorry.", m);
			disp(str)
	end
end
%=======================================================================================================================

%=======This function ensures that the inpts are of allowable type/value/etc.===========================================
function check_inputs(v, w, m, k)
	validateattributes(v, {'numeric'}, {'scalar', 'real'})			%Real number
	validateattributes(w, {'numeric'}, {'scalar', 'real'})			%Real number
	validateattributes(m, {'numeric'}, {'scalar', 'integer', 'positive'})	%Real number
	validateattributes(k, {'numeric'}, {'vector', 'real'})			%Vector of real numbers
end
%=======================================================================================================================


