%	plot_state_on_lattice
%This function plots an eignestate of a lattice Hamiltonian on the lattice itself. We need the coordinates of each 
%lattice point as well as the eigenstate in question as a vector. 
%
%Inputs:
%	X	-	a NxN square matrix of real values, where N is the number of sites, listing x-coordinate values
%	Y	-	a NxN square matrix of real values listing y-coordinate values
%	evect	-	a 1xN or Nx1 complex vector representing the eigenstate overlap with each lattice site
%	npts	-	number of data visualization points in each direction (square); must be positive integer
%	fig	-	a figure handle (CAUTION!! This figure will be cleared!)
%
%Outputs:
%	data_out	-	the data plotted on fig, an npts x npts array of positive real values
%
%Example Usage:
%	data_out = plot_state_on_lattice(X, Y, evect, npts, fig)
%
function data_out = plot_state_on_lattice(X, Y, evect, npts, fig)
	check_inputs(X, Y, evect, npts, fig)					%Make sure inputs are of allowable types
	X = diag(X); X = full(X); Y = diag(Y); Y = full(Y);			%Fix formatting
	num_sites = length(evect); 						%Number of sites
	[x_coords, y_coords] = set_up_coordinates(X, Y, npts);			%Get array of coordinates for plotting
	ax = set_up_plot(fig, x_coords, y_coords);				%Clear figure and plot lattice
	data_out = zeros(npts, npts);						%Will hold data for plotting
	width = 0.25; 								%Hardcoded gaussian width parameter
	for index = 1:num_sites							%Loop over all sites
		strength = realsqrt(evect(index)' * evect(index)); 		%Magnitude of amplitude at lattice site
		data_out = data_out + gaussian(x_coords, y_coords,...		%Add one gaussian to data
				X(index), Y(index), strength, width, npts); 
	end
	put_data_on_plot(ax, x_coords, y_coords, data_out)			%Draw the data
end

%=======check_inputs====================================================================================================
%Makes sure that input values are of allowable types and fall within acceptable ranges
function check_inputs(X, Y, evect, npts, fig)
	validateattributes(X, {'numeric'}, {'2d', 'square', 'nonempty', 'real'})%A non-empty square matrix of reals
	validateattributes(Y, {'numeric'}, {'2d', 'square', 'nonempty', 'real'})%A non-empty square matrix of reals
	validateattributes(evect, {'numeric'}, {'vector', 'nonempty'})		%A non-empty vector of numbers
	assert(size(X, 1) == length(evect) && size(Y, 1) == length(evect))	%All lengths must be compatible
	offDiagMask = ~eye(size(X));						%Make sure that X and Y are diagonal
	assert(~any(abs(X(offDiagMask)) > 1e-12))
	assert(~any(abs(Y(offDiagMask)) > 1e-12))
	validateattributes(npts, {'numeric'}, {'scalar', 'positive', 'integer'})%Positive integer
	%The two-step verification is necessary to avoid a weird corner case.
	%See https://www.mathworks.com/matlabcentral/answers/300880-what-is-best-practice-to-determine-if-input-is-a-figure-or-axes-handle.
	assert(ishghandle(fig, 'figure') && isa(fig, 'matlab.ui.Figure')) 
end
%=======================================================================================================================

%=======set_up_plot=====================================================================================================
%Creates the axes on fig and does axis label and title formatting. Returns axes object. 
function ax = set_up_plot(fig, x_coords, y_coords)
	clf(fig); 								%Clear the figure of EVERYTHING!
	ax = axes(fig); 							%Create axes object on fig
	axis(ax, 'xy');
	padding_factor = 1.0;
	xlimmin = padding_factor * min(x_coords); xlimmax = padding_factor * max(x_coords);
	ylimmin = padding_factor * min(y_coords); ylimmax = padding_factor * max(y_coords);
	%k_x axis setup								%Standard settings for this kind of plot
		xlim(ax, [xlimmin, xlimmax]);
		%set(ax, 'XTick', [-pi, -pi/2, 0, pi/2, pi]);
		%set(ax, 'XTickLabels', ["$-\frac{\pi}{a}$", "$-\frac{\pi}{2a}$", "$0$", "$\frac{\pi}{2a}$", "$\frac{\pi}{a}$"]);
		%xlabel(ax, "$k_x$", 'interpreter', 'latex', 'rotation', 0);
		xlabel(ax, "$x$", 'interpreter', 'latex', 'rotation', 0);
	%k_y axis setup
		ylim(ax, [ylimmin, ylimmax]);
		%set(ax, 'YTick', [-pi, -pi/2, 0, pi/2, pi]);
		%set(ax, 'YTickLabels', ["$-\frac{\pi}{a}$", "$-\frac{\pi}{2a}$", "$0$", "$\frac{\pi}{2a}$", "$\frac{\pi}{a}$"]);
		%ylabel(ax, "$k_y$", 'interpreter', 'latex', 'rotation', 0);
		ylabel(ax, "$y$", 'interpreter', 'latex', 'rotation', 0);
	ax.TickLabelInterpreter = 'latex';
	ax.FontSize = 12;
	ax.XRuler.TickLength = [0, 0]; 
	ax.YRuler.TickLength = [0, 0]; 
	str = 'Amplitude of Approximate Eigenstate at Lattice Sites';
	title(ax, str, 'interpreter', 'latex', 'FontSize', 16);
	pbaspect(ax, [1,1,1])
	hold(ax, 'on')								%Hold for later data
end
%=======================================================================================================================

%=======set_up_coordinates==============================================================================================
%Create the lists of coordinates for plotting
function [x_coords, y_coords] = set_up_coordinates(X, Y, npts);			%Get array of coordinates for plotting
	padding_factor = 1.1;
	xmin = padding_factor * min(X(:)); xmax = padding_factor * max(X(:));	%Extreme values of lattice coords
	ymin = padding_factor * min(Y(:)); ymax = padding_factor * max(Y(:));
	x_coords = linspace(xmin, xmax, npts); 					%Make the lists
	y_coords = linspace(ymin, ymax, npts); 
end
%=======================================================================================================================

%=======gaussian========================================================================================================
%Computes the values of a gaussian function with given strength, width, and center
function out = gaussian(x_coords, y_coords, x_center, y_center,...
							strength, width, npts)
	out = zeros(npts, npts);						%Initialize output
	for index = 1:npts
		for jndex = 1:npts
			sq_dist = (x_coords(index) - x_center)^2 + ...		%Distance from lattice site
				(y_coords(jndex) - y_center)^2; 
			st_dev = 2 * width^2; 					%Standard deviation
			out(jndex, index) = strength * exp(-sq_dist) / st_dev;	%Formula for gaussian
		end
	end
end
%=======================================================================================================================

%=======put_data_on_plot================================================================================================
%Draws the actual data
function put_data_on_plot(ax, x_coords, y_coords, data_out)					
	heatmap = imagesc(ax, x_coords, y_coords, data_out); 			%Create heatmap
%	colormap(ax, flipud(cmocean('deep')));		%A good perceptually uniform colormap, but reversed!!
	colormap(ax, cmocean('deep'));
%	clim(ax, [-2, 1])
%	clim(ax, color_lim)
	cb = colorbar(ax);
%	cb.Ticks = linspace(color_lim(1), color_lim(2), 5);
	hold(ax, 'off')								%Stop holding for more data
end
%=======================================================================================================================
