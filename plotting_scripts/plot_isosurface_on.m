%=======plot_isosurface_on==============================================================================================
%
%This function plots isosurfaces of some 3d data at specified iso values on a given figure object.
%WARNING--the figure handle input to this function will be fully cleared. 
%
%Inputs:
%	fig   		- 	figure handle (e.g., figure(1)) THIS IS FULLY CLEARED
%   	data  		- 	numeric 3D array of reals
%   	iso_vals 	- 	vector of real scalar isovalues
%	view_angle 	-	a 1x3 vector pointing in the direction of the camera
%
%Outputs:
%	patch_obj 	-	a patch object approximating the isosurface
%
%Example Usage:
%	patch_obj = plot_isosurface_on(fig, data, iso_vals, [1, 1, 1])
function patch_obj = plot_isosurface_on(fig, data, iso_vals, view_angle)
	check_inputs(fig, data, iso_vals, view_angle)			%Make sure inputs are of okay types, sizes, etc. 
	[ax, colors] = axis_setup(fig, view_angle, iso_vals);		%Set up the axis object inside fig
	patch_obj = add_data_to_plot(ax, data, iso_vals, colors);	%Draw isosurfaces on plot
end
%=======================================================================================================================

%Makes sure that inputs are of correct types and sizes
function check_inputs(fig, data, iso_vals, view_angle)
	%The two-step verification is necessary to avoid a weird corner case.
	%See https://www.mathworks.com/matlabcentral/answers/300880-what-is-best-practice-to-determine-if-input-is-a-figure-or-axes-handle.
	assert(ishghandle(fig, 'figure') && isa(fig, 'matlab.ui.Figure')) 
	validateattributes(data, {'numeric'}, {'real', '3d'})		%A 3D array of real values
	validateattributes(iso_vals, {'numeric'}, {'vector', 'real'})	%A real vector 
	validateattributes(view_angle, {'numeric'}, {'vector', 'size', [1, 3], 'real'}) %Something like [1, 1, 1]
end
%=======================================================================================================================

%=======axis_setup======================================================================================================
%Sets up an axis object inside the figure and returns the axis object, as well as colors for use in plotting
function [ax, colors] = axis_setup(fig, view_angle, iso_vals, npts_array)

	cmap = flipud(cmocean('deep'));	%A good perceptually uniform colormap
	vmin = min(iso_vals); 		%Normalize isovalues to [0,1] using their min/max 
	vmax = max(iso_vals);
	t = (iso_vals - vmin) ./ max(eps, (vmax - vmin));   % normalized [0,1]
	M = size(cmap,1); 		%Map normalized positions to colormap indices (1..M) and interpolate
	idx = 1 + (M-1)*t;              %fractional indices in 1..M
	r = interp1(1:M, cmap(:,1), idx);
	g = interp1(1:M, cmap(:,2), idx);
	b = interp1(1:M, cmap(:,3), idx);
	colors = [r(:), g(:), b(:)];   %N×3 RGB colors for each iso_val

    	clf(fig);						%Clear the figure
	ax = axes('Parent', fig);				%Create axes
	hold(ax, 'on');						%Plot multiple things without erasing
	%axis(ax, 'vis3d');					%It's a 3d plot
	box(ax, 'on');				
	grid(ax, 'on');
	daspect(ax, [1 1 1]);					%Aspect ratio
	camlight(ax,'headlight');				%Good lighting options
	lighting(ax,'gouraud');
	material(ax,'dull');
	%k_x axis setup
		xlim(ax, [-pi, pi]);
		set(ax, 'XTick', [-pi, -pi/2, 0, pi/2, pi]);
		set(ax, 'XTickLabels', ["$-\frac{\pi}{a}$", "$-\frac{\pi}{2a}$", "$0$", "$\frac{\pi}{2a}$", "$\frac{\pi}{a}$"]);
		xlabel(ax, "$k_x$", 'interpreter', 'latex', 'rotation', 0);
	%k_y axis setup
		ylim(ax, [-pi, pi]);
		set(ax, 'YTick', [-pi, -pi/2, 0, pi/2, pi]);
		set(ax, 'YTickLabels', ["$-\frac{\pi}{a}$", "$-\frac{\pi}{2a}$", "$0$", "$\frac{\pi}{2a}$", "$\frac{\pi}{a}$"]);
		ylabel(ax, "$k_y$", 'interpreter', 'latex', 'rotation', 0);
	%E axis setup
		zlim(ax, [-4, 4])
		set(ax, 'ZTick', [-4, -2, 0, 2, 4]);
		set(ax, 'ZTickLabels', ["-4", "-2", "0", "2", "4"]);
		zlabel(ax, "$E$", 'interpreter', 'latex', 'rotation', 0);
	ax.TickLabelInterpreter = 'latex';
	ax.FontSize = 12;
	ax.XRuler.TickLength = [0, 0]; 
	ax.YRuler.TickLength = [0, 0]; 
	str = '$\mu^Q_{H, T_1, T_2, T_3}$';
	title(ax, str, 'interpreter', 'latex', 'FontSize', 16);
	view(ax, view_angle);					%View angle

	% Create a mapped colormap and colorbar
	colormap(ax, cmap);
	c = colorbar(ax);
	caxis(ax, [vmin vmax]);                       % map data range to colormap
	% Place ticks at your isovalues (if too many, pick a subset)
	c.Ticks = iso_vals;                           % requires iso_vals within [vmin,vmax]
	c.TickLabels = arrayfun(@num2str, iso_vals, 'UniformOutput', false);
end
%=======================================================================================================================

%=======add_data_to_plot================================================================================================
%Draws the actual isosurface on the plot and returns the patch object
function patch_objs = add_data_to_plot(ax, data, iso_vals, colors)
	patch_objs = { }; 					%Initially empty
	x_vals = linspace(-pi, pi, size(data, 1));		%Create arrays of x, y, and E values 
	y_vals = linspace(-pi, pi, size(data, 2)); 
	E_vals = linspace(-4, 4, size(data, 3)); 
	for index = 1:length(iso_vals)
		fv = isosurface(x_vals, y_vals, E_vals, data, iso_vals(index));	%Get isosurface 
		patch_objs{end + 1} = patch(ax, fv,...		%Create patch object and add to array
					'FaceColor', colors(index, :),...
					'EdgeColor', 'none',...
					'FaceAlpha', 1.0);	
	end
end
%=======================================================================================================================
