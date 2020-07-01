function [] = function_3Dfigure_withpatch(data_figure,parameters_figure)

Fig_ = figure; % Create figure
Fig_.Name= parameters_figure.figurename;
Fig_.Color='white'; % Background colour
set(Fig_, 'Position', parameters_figure.figureposition)
axes_ = axes('Parent',Fig_); % Create axes
hold(axes_,'on');
t_=title (' ','FontName',parameters_figure.fontname,'FontSize',16);  % Title
t_.String= parameters_figure.title;
patch('Faces',data_figure.face,'Vertices',data_figure.vertices,'FaceVertexCData',data_figure.FaceVertexCData,'FaceColor',data_figure.FaceColor,'EdgeColor',data_figure.EdgeColor)
axis equal; axis tight;
xlabel(parameters_figure.xlabel); % - Axis label
ylabel(parameters_figure.ylabel);
zlabel(parameters_figure.zlabel);
grid(axes_,parameters_figure.grid); % Display grid
set(axes_,'XMinorGrid',parameters_figure.minorgrid,'YMinorGrid',parameters_figure.minorgrid,'ZMinorGrid',parameters_figure.minorgrid); % Display grid for minor thicks also
set(axes_,'FontName',parameters_figure.fontname,'FontSize',14); % Fontname and fontsize
%axis vis3d tight % view;
camlight left
%daspect([1,1,1])
view(3)
colormap(parameters_figure.colormap); % Color map
lighting(parameters_figure.lighting) % Light
colorbar('peer',axes_); % Create colorbar
hold(axes_,'off');
% Save figures
filename = parameters_figure.filename;
fullpath = parameters_figure.fullpath;
savefig(Fig_,[fullpath filename]) % .fig
saveas(Fig_,[fullpath filename],'png') % .png
close(Fig_)

end

