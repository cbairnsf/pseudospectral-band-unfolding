%	maybe_create_parpool 
%This function creates a parpool only if one does not already exist and sets parfor opts.
%The parpool size is one less than the number of CPU cores available, so that one is free for other tasks. 
%The 'RangePartitionMethod' is set to 'fixed' and the 'SubrangeSize' is set to 2; those choices have been working. 
%(It might be worthwhile to play around with those options to try to optimize performance.)
%Warnings are disabled on the parallel workers. 
%NB the warning 'MATLAB:datetime:NonstandardSystemTimeZoneFixed' is supressed as well (permanently!)
%
%	Warnings:
%		This requires a UNIX environment to work, since it makes an 'nproc' call to the system. 
%		That will fail on Windows. 
%		NB - warnings will be disabled for all parallel workers if a new parpool is created. 
%
%	Outputs: 
%		parpool_handle -- a handle to a parallel pool 
%		opts -- the options object (output of parforOptions)
%
%	Example Usage:
%		[parpool_handle, opts] = maybe_create_parpool();
function [parpool_handle, opts] = maybe_create_parpool()
	warning('off','MATLAB:datetime:NonstandardSystemTimeZoneFixed')	%This warning is annoying
	parpool_handle = gcp('nocreate');			%Check for existing parpool
	if isempty(parpool_handle)	
		[~, numprocs] = system('nproc');		%Doing this in two lines avoids shell output
		workers = str2num(numprocs) - 1;		%The number of parallel workers--leave one free
		parpool_handle = parpool("Processes", workers);
		pctRunOnAll warning off;			%Disable warnings for parallel workers
	end
	opts = parforOptions(parpool_handle,'RangePartitionMethod','fixed','SubrangeSize',2);
end
