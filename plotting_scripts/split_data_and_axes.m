%	split_data_and_axes.m
%This function splits the data and axes of a figure into separate figures so that they can be saved separately. 
%
%Inputs:
%	fig0	-	a figure handle object 
%	fig1	-	a figure handle object (NB this will be cleared!)
%	fig2	-	a figure handle object (NB this will be cleared!)
%
%Outputs:
%	fig0 is unchanged
%	fig1 contains axes but no data (to be saved as an svg)
%	fig2 contains data but no axes (to be saved somehow)
%
%Example Usage:
%	split_data_and_axes(fig0, fig1, fig2);
%
function split_data_and_axes(fig0, fig1, fig2)
	clf(fig1)				%Clear the figures that will hold the copies
	clf(fig2)
	copyobj(allchild(fig0), fig1);		%Copy image twice and create axis hanles
	ax1 = fig1.Children(end);
	copyobj(allchild(fig0), fig2);
	ax2 = fig2.Children(end);

	%Create the image with only axes
	lines = findobj(ax1, 'Type','line');	%Make all data objects invisible
	patches = findobj(ax1, 'Type','patch');
	surfs = findobj(ax1, 'Type','surface');
	images = findobj(ax1, 'Type','image');
	set([lines; patches; surfs; images], 'Visible', 'off');
	%delete(fig1.Children(1)); 		%Delete the colorbar
	

	%Create the image with only data
	set(fig2, 'Color', 'none');		%Transparent background
	set(ax2, 'Color', 'none');

	if isprop(ax2, 'Title'),		%Hide title and labels
		ax2.Title.String = '';
	end
	ax2.XLabel.String = '';
	ax2.YLabel.String = '';
	if isprop(ax2, 'ZLabel')
		ax2.ZLabel.String = '';
	end

	ax2.Box = 'off';			%Hide axes completely
	ax2.XColor = 'none';
	ax2.YColor = 'none';
	if isprop(ax2,'ZColor')
		ax2.ZColor = 'none';
	end
	ax2.XTick = [];
	ax2.YTick = [];
	if isprop(ax2, 'ZTick')
		ax2.ZTick = [];
	end

	grid(ax2, 'off');			%Hide grids
	if isprop(ax2, 'MinorGridLineStyle')
		ax2.MinorGridLineStyle = 'none';
	end
	%delete(fig2.Children(1)); 		%Delete the colorbar
end
