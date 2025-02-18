function [] = Function_plotcorrelationheatmap(CorrMatrix,Parnames,metric,plotper,pformat)

if strcmp(metric,'Kendall') % Need diverging colormap
    cmap = crameri('vik'); % CALL THIRD PARTY COLOR MAP
    % Crameri, F. (2018). Scientific colour-maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
    % Crameri, F. (2018), Geodynamic diagnostics, scientific visualisation and StagLab 3.0, Geosci. Model Dev., 11, 2541-2562, doi:10.5194/gmd-11-2541-2018.
    % For more on choosing effective and accurate colormaps for science, be sure to enjoy this fine beach reading:
    % Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco. 2016. True colors of oceanography: Guidelines for effective and accurate colormap selection. Oceanography 29(3):9-13, http://dx.doi.org/10.5670/oceanog.2016.66.
end

Fig = figure;
Fig.Name= ['Correlation per ' plotper];
Fig.Color='white'; % Background colour
%[t,s] = title('Correlation heatmap',[plotper ', ' metric]);
%t.FontSize = pformat.size.sgtitle_main; s.FontSize = pformat.size.sgtitle_sub;
%s.FontAngle = 'italic';
h=heatmap(Parnames,Parnames,CorrMatrix,'Colormap',cmap,'ColorLimits',[-1 1],'GridVisible','off');
if ~pformat.showcellvalue 
    h.CellLabelColor = 'none';
end
h.Title = {'Correlation heatmap',['Per ' plotper ', ' metric]};
h.FontName = pformat.fontname;
h.FontSize = pformat.size.axes;

end