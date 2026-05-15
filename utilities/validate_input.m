%	validate_input
%This function tests to make sure a variable is a numeric scalar and obeys certain restrictions.
%
%	Input (repeated any number of times): 
%		var -- some variable
%		srt -- a string describing the requirements
%
%	Options for str:
%		"numeric"
%		"real"
%		"non-negative real"
%		"positive real"
%		"integer"
%		"non-negative integer"
%		"positive integer"
%
%	Output:
%		truth_value -- this is 1 if each variable passes the test and 0 otherwise
%
%	Example usage:
%		truth_value = validate_input(var, "non-negative integer")
function truth_value = validate_input(var, str)
	truth_value = 0;		%Default to failing
	if ~isstring(str)		%First make sure this function is being used correctly
		error("second input must be a string")	%Bail--the function has been called wrongly	
	elseif ~(isnumeric(var) && isscalar(var))
		%do nothing--var has failed
	else	%only get to here if var is an array of numeric type
		switch str		%Test the input variable based on the requirement string
			case "numeric"
				truth_value = 1; %We already passed!
			case "real"
				if all(isreal(var(:)))	%The var(:) turns it into a vector so all gives a scalar
					truth_value = 1;
				end
			case "non-negative real"
				if all(isreal(var(:))) && all(var(:) >= 0)
					truth_value = 1;
				end
			case "positive real"
				if all(isreal(var(:))) && all(var(:) > 0)
					truth_value = 1;
				end
			case "integer"
				if all(var(:) == floor(var(:)))
					truth_value = 1;
				end
			case "non-negative integer"
				if all(var(:) == floor(var(:))) && all(var(:) >= 0)
					truth_value = 1;
				end
			case "positive integer"
				if all(var(:) == floor(var)) && all(var(:) > 0)
					truth_value = 1;
				end
			otherwise 
				error("invalid string was passed")	%If the second input is not from list
		end
end
