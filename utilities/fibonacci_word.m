%	fibonacci_word
%
%This function produces the word corresponding to the nth Fibonacci number. The algorithm is to start with B and then
%apply the substitution rules
%	B --> A
%	A --> AB
%n times. We start with the word "AAB" (Terry's modification). The resulting word will *not* have length F_n. Some examples:
%
%	n 	C_n 		F_n
%	0 	B 		1
%	1 	A 		1
%	2 	AB 		2
%	3 	ABA 		3
%	4 	ABAAB 		5
%	5 	ABAABABA 	8
%	6 	ABAABABAABAAB 	13
%
%Inputs:
%	n 	--	the number of substitutions
%	A	-- 	the letter to use for A (single character)
%	B	-- 	the letter to use for B (single character, should ideally be different!) 
%	repl	-- 	what is A replaced with in each iteration (a character array)
%
%Outputs:
%	C_n	-- 	the nth Fibonacci word using the letters A and B 
%
%Example Usage:
%	word = fibonacci_word(64, v, w, 'vwww'); 
%
%repl options:
%	'AB'	-- 	the ``golden'' option (fibonacci)
%	'AAB'	--	the ``silver'' option
%	'ABBB'	--	non quasi-crystal considered in reference
%
%Reference:
%	Jagannathan, A. (2021).
%	The Fibonacci quasicrystal: Case study of hidden dimensions and multifractality.
%	Reviews of Modern Physics, 93(4), 045001.
%	https://doi.org/10.1103/RevModPhys.93.045001
%
function word = fibonacci_word(n, A, B, repl)
	check_inputs(n, A, B, repl);		%Make sure inputs are allowable
	word = strcat(A, A, B);			%Initial state	(This is the Terry modification.) 
	for index = 1:n				%The recursive loop
		word_modified = ''; 		%Empty word
		for jndex = 1:length(word)
			if word(jndex) == A
				word_modified = [word_modified, repl];
			else
				word_modified = [word_modified, A];  
			end
		end		
		word = word_modified; 
	end
end

%This function makes sure the inputs are valid
function check_inputs(n, A, B, repl)
	validateattributes(n, {'numeric'}, {'scalar', 'nonnegative', 'integer'})	%A positive integer
	validateattributes(A, {'char'}, {'scalartext', 'nonempty', 'numel' 1})		%A one entry character array 
	validateattributes(B, {'char'}, {'scalartext', 'nonempty', 'numel' 1})		%A one entry character array 
	validateattributes(repl, {'char'}, {'scalartext', 'nonempty'})			%A one entry character array a
	for index = 1:length(repl)
		if (repl(index) ~= A) && (repl(index) ~= B) 				%Make sure no weird letters in repl
			error("Extraneous character in replacement word.")
		end
	end
end
