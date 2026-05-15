%	update_progress_bar
%This functions adds one pipe ('|') to the repl. It is meant to be called repeatedly in a loop to act as a progress bar. 
%The progress bar does not work properly in a terminal, hence the wrapping in the if statement in the example below. 
%Visually, you should see this:
%	..............................
%	||||||||||||||
%
%Making this a separate function is a bit silly, but it's now "self-documenting code", the \b|\n stuff is hard to parse. 
%
%	Example usage:
%		prog_bar_exists = maybe_create_progress_bar(num_iterations);
%		for index = 1:N
%			<loop code here>
%			if prog_bar_exists
%				update_progress_bar()
%			end
%		end
function update_progress_bar()
	fprintf('\b|\n');	%\b deletes \n and returns us to the previous line, then we draw '|' and return
end
