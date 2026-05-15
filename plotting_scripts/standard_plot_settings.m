%	standard_plot_settings
%This function applies standard settings to a heatmap to save typing.
%
%	Inputs:
%		ax -- an axes handle;
%		color_lim -- a 1x2 vector of positive real numbers to be the colorbar limits
%		heatmap_index -- which child of ax is the heatmap? ax.Children(heatmap_index) should be the heatmap.
%		reverse_flag -- "r" to reverse the colormap
%
%	fig should have been created via something like fig = imagesc(data_out);
%	if there are more 'Children' under fig then there could maybe be problem?
%
%	Example usage:
%		standard_plot_settings(ax, color_lim, 1, 'r');
function standard_plot_settings(ax, color_lim, heatmap_index, varargin)
	%Hardcoded
	%fs = 9; 					%The font size
	fs = 20;
	ncbt = 5;					%Number of colorbar tick marks
	if nargin >= 4 & strcmp(varargin{1}, "r")
		reverse_flag = 1;
	else
		reverse_flag = 0;
	end
	axis(ax, "xy");					%Larger values of y on top
	set(ax.Parent, 'InvertHardcopy', 'off',...	%Don't swap white/black ever
		'Color', [1, 1, 1])			%Pure white
	set(ax, 'XColor', [0, 0, 0],...			%Pure black
		'YColor', [0, 0, 0],...		
		'ZColor', [0, 0, 0],...	
		'GridColor', [200, 200, 200]/256,...	%Grey
		'GridAlpha', 1)			%Fully opaque

	heatmap = ax.Children(heatmap_index);	%Handle to image
	npts = size(heatmap.CData);		%Get size of data
	xvals = linspace(-pi, pi, npts(2));	%Values along horizontal axis
	yvals = linspace(-4, 4, npts(1));	%Values along vertical axis
	%yvals = linspace(-2.0, 2.0, npts(1));	%Values along vertical axis
	heatmap.XData = xvals;			%Set x,y values correctly
	heatmap.YData = yvals;
	%Slice through 3d version:
	%	xlim(ax, [-pi, pi]);
	%	set(ax, 'XTick', [-pi, -pi/2, 0, pi/2, pi]);
	%	set(ax, 'XTickLabels', ["$-\frac{\pi}{a}$", "$-\frac{\pi}{2a}$", "$0$", "$\frac{\pi}{2a}$", "$\frac{\pi}{a}$"]);
	%	xlabel(ax, "$k_x$", 'interpreter', 'latex', 'rotation', 0);
	%	ylim(ax, [-pi, pi]);
	%	set(ax, 'YTick', [-pi, -pi/2, 0, pi/2, pi]);
	%	set(ax, 'YTickLabels', ["$-\frac{\pi}{a}$", "$-\frac{\pi}{2a}$", "$0$", "$\frac{\pi}{2a}$", "$\frac{\pi}{a}$"]);
	%	ylabel(ax, "$k_y$", 'interpreter', 'latex', 'rotation', 0);
	%2d version:
		xlim(ax, [-pi, pi]);
		set(ax, 'XTick', [-pi, -pi/2, 0, pi/2, pi]);
		set(ax, 'XTickLabels', ["$-\frac{\pi}{a}$", "$-\frac{\pi}{2a}$", "$0$", "$\frac{\pi}{2a}$", "$\frac{\pi}{a}$"]);
		xlabel(ax, "$k$", 'interpreter', 'latex', 'rotation', 0);
		ylim(ax, [-4, 4]);
		%ylim(ax, [-2, 2]);
		set(ax, 'YTick', [-4, -2, 0, 2, 4]);
		set(ax, 'YTickLabels', ["-4", "-2", "0", "2", "4"]);
		%set(ax, 'YTick', [-2, -1, 0, 1, 2]);
		%set(ax, 'YTickLabels', ["-2", "-1", "0", "1", "2"]);
		ylabel(ax, "$E$", 'interpreter', 'latex', 'rotation', 0);
	%title(ax, "$-\log\mu^Q_{H, T} (E, k)$", 'interpreter', 'latex', 'FontSize', 20);	%Set the title
	%title(ax, "$\mu^Q_{H, T} (E, k)$", 'interpreter', 'latex', 'FontSize', fs); 
	%title(ax, "percent difference $\frac{OBC - PBC}{OBC}$", 'interpreter', 'latex', 'FontSize', fs); 
	title(ax, "$\mu^M_{H, T} (E, k)$", 'interpreter', 'latex', 'FontSize', fs, 'Color', 'k');
	ax.TickLabelInterpreter = 'latex';
	ax.FontSize = fs;
	ax.XRuler.TickLength = [0, 0]; 
	ax.YRuler.TickLength = [0, 0]; 
%	title(ax, "$-\log\mu^Q_{H, T, X} (E, k, x)$", 'interpreter', 'latex', 'FontSize', 20);	%Set the title
	ax.TickLabelInterpreter = 'latex';
	ax.FontSize = fs;
	ax.XRuler.TickLength = [0, 0]; 
	ax.YRuler.TickLength = [0, 0]; 
%	title(ax, "$-\log\mu^Q_{H, T, X} (E, k, x)$", 'interpreter', 'latex', 'FontSize', 20);	%Set the title
	if reverse_flag
		colormap(ax, flipud(cmocean('deep')));		%A good perceptually uniform colormap, but reversed!!
	else
		colormap(ax, cmocean('deep'));			%A good perceptually uniform colormap
	end
	clim(ax, color_lim)
	cb = colorbar(ax);
	cb.Limits = color_lim;
	cb.Ticks = linspace(color_lim(1), color_lim(2), ncbt);
	cb.Ruler.TickLabelFormat = '%0.3f';
	cb.Color = 'k';
end
