%	maybe_create_progress_bar
%This function creates a progress bar for use during a loop, but only in the desktop environment, not the cli. 
%(The progress bar does not work properly in a terminal.) 
%Visually, you should see this:
%	..............................
%	||||||||||||||
%
%	Input:
%		num_iterations -- a non-negative integer
%	Output:
%		a row of num_interations dots are printed to the repl
%		prog_bar_exists -- 1 if a progress bar was created (so that it can be updated later) or a 0 otherwise
%
%	Example usage:
%		prog_bar_exists = maybe_create_progress_bar(num_iterations);
%		for index = 1:N
%			<loop code here>
%			if prog_bar_exists
%				update_progress_bar()
%			end
%		end
%
function prog_bar_exists = maybe_create_progress_bar(num_iterations)
	prog_bar_exists = 0;						%Default to no bar
	if ~validate_input(num_iterations, "non-negative integer")	%Test for valid input
		error("the number of iterations must be a non-negative integer")
	elseif usejava('desktop')					%Test for desktop environment
		fprintf(['\n' repmat('.',1,num_iterations) '\n\n']);	%Display a progress bar then two blank lines
		prog_bar_exists = 1;					%We made a progress bar!
	end
end
