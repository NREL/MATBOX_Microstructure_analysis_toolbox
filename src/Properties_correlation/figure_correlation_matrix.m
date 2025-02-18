function figure_correlation_matrix(parms, volume_choice, group_choice, plotper, pformat)

markerlist = {'o','square','diamond','^','v','<','pentagram','hexagram','+','*','.','x','_','|'};
markerlist = [markerlist markerlist markerlist markerlist markerlist];
markerlist = [markerlist markerlist markerlist markerlist markerlist];
colorlist = [colororder; rand(100,3)];

[number_parameter,number_volume] = size(parms);
number_volume = number_volume-2;
[number_group,~] = size(group_choice);

plotdist = false; % Beta

% Clean vars
tmp = parms(:,3:end);
parmsvar = NaN(number_parameter,number_volume);
for kpar=1:1:number_parameter
    for kvol=1:1:number_volume
        v = tmp(kpar,kvol).Variables;
        if iscell(v)
            v = cell2mat(v);
            if ~isempty(v)
                parmsvar(kpar,kvol) = v;
            end
        else
            parmsvar(kpar,kvol) = v;
        end
    end
end

% Group id
if strcmp(plotper,'group')
    group_id = cell2mat(group_choice(:,1));
    volume_group = cell2mat(volume_choice(:,3));
end

Fig = figure; % Create figure
Fig.Name= ['Correlation per ' plotper];
Fig.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig,'position',[scrsz(1) scrsz(2) 0.9*scrsz(3) 0.8*scrsz(4)]);

t = tiledlayout(number_parameter,number_parameter,'TileSpacing','Compact','Padding','Compact');
title(t,'Correlation matrix','FontWeight','normal','Fontname',pformat.fontname,'Fontsize',pformat.size.sgtitle_main)
t.Subtitle.String = ['Per ' plotper];
t.Subtitle.FontAngle = 'italic';
t.Subtitle.FontName = pformat.fontname;
t.Subtitle.FontSize = pformat.size.sgtitle_sub;

str_eval = '';
for row=1:1:number_parameter
    for col=1:1:number_parameter
        nexttile
        if col==row

            % Set axis label
            if (col==1 && row==1 ) || (col==number_parameter && row==number_parameter)
                if col==1 && row==1
                    par_name = char(parms.Parameter(1));
                    unit_name = char(parms.Unit(1));
                else
                    par_name = char(parms.Parameter(number_parameter));
                    unit_name = char(parms.Unit(number_parameter));
                end
                if isempty(unit_name)
                    str_axis = par_name;
                else
                    if strcmp(unit_name,'um') || strcmp(unit_name,'micrometer') || strcmp(unit_name,'Micrometer') || strcmp(unit_name,'micrometers') || strcmp(unit_name,'Micrometers')
                        unit_name = '(\mum)';
                    elseif strcmp(unit_name,'um-1') || strcmp(unit_name,'micrometer-1') || strcmp(unit_name,'Micrometer-1') || strcmp(unit_name,'micrometers-1') || strcmp(unit_name,'Micrometers-1')
                        unit_name = '(\mum^{-1})';
                    elseif strcmp(unit_name,'um-2') || strcmp(unit_name,'micrometer-2') || strcmp(unit_name,'Micrometer-2') || strcmp(unit_name,'micrometers-2') || strcmp(unit_name,'Micrometers-2')
                        unit_name = '(\mum^{-2})';                        

                    elseif strcmp(unit_name,'nm') || strcmp(unit_name,'nanometer') || strcmp(unit_name,'Nanometer') || strcmp(unit_name,'nanometers') || strcmp(unit_name,'Nanometers')
                        unit_name = '(nm)';
                    elseif strcmp(unit_name,'nm-1') || strcmp(unit_name,'nanometer-1') || strcmp(unit_name,'Nanometer-1') || strcmp(unit_name,'nanometers-1') || strcmp(unit_name,'Nanometers-1')
                        unit_name = '(nm^{-1})';
                    elseif strcmp(unit_name,'nm-2') || strcmp(unit_name,'nanometer-2') || strcmp(unit_name,'Nanometer-2') || strcmp(unit_name,'nanometers-2') || strcmp(unit_name,'Nanometers-2')
                        unit_name = '(nm^{-2})';   

                    else
                        unit_name = ['(' unit_name ')'];
                    end
                    str_axis = [par_name ' ' unit_name];
                end
                if col==1 && row==1
                    h(1) = ylabel(str_axis);
                else
                    h(1) = xlabel(str_axis);
                end
                set(gca, 'Xcolor', 'w', 'Ycolor', 'w');
                set(h, 'Color', 'k');
                set(gca,'Fontname',pformat.fontname,'Fontsize',pformat.size.axes)      
            else
                set(gca,'Visible','off')
            end

        else

            if strcmp(plotper,'volume') % Plot per volume
                hold on;
                for kvol = 1:1:number_volume
                    x = parmsvar(col,kvol);
                    y = parmsvar(row,kvol);
                    if ~isnan(x) && ~isnan(y)
                        if sum(isnan(parmsvar(col,:))) || sum(isnan(parmsvar(row,:)))
                            plot(x,y,"LineStyle","none","Marker",markerlist(kvol),"MarkerSize",pformat.markersize,"LineWidth",pformat.linewidth,"Color",colorlist(kvol,:));
                        else
                            hvol(kvol).plot = plot(x,y,"LineStyle","none","Marker",markerlist(kvol),"MarkerSize",pformat.markersize,"LineWidth",pformat.linewidth,"Color",colorlist(kvol,:),'DisplayName',char(parms.Properties.VariableNames(kvol+2)));
                        end
                        if pformat.text
                            text(x,y,char(parms.Properties.VariableNames(kvol+2)),'Fontname',pformat.fontname,'Fontsize',pformat.textsize);
                        end
                    end  
                end
                hold off;
                
            elseif strcmp(plotper,'group') % Plot per group
                if ~plotdist
                    hold on;
                else
                    xvalues = []; yvalues = []; grpvalues = {}; scol = [];
                end

                for kgr = 1:1:number_group
                    idg = group_id(kgr);
                    idv = find(volume_group==idg);
                    if ~isempty(idv)
                        xy = [parmsvar(col,idv)' parmsvar(row,idv)'];
                        if pformat.text
                            [n,~]=size(xy);
                            for kv=1:1:n
                                text(xy(kv,1),xy(kv,2),char(parms.Properties.VariableNames(idv(kv)+2)),'Fontname',pformat.fontname,'Fontsize',pformat.textsize);
                            end
                        end
                        xy = sortrows(xy,1); % Useful is LineStyle is not none
                        if ~plotdist
                            if cell2mat(group_choice(kgr,8))
                                hvol(kgr).plot = plot(xy(:,1),xy(:,2),"LineStyle",char(group_choice(kgr,4)),"Marker",char(group_choice(kgr,6)),"MarkerSize",cell2mat(group_choice(kgr,7)),"MarkerFaceColor",str2num(cell2mat(group_choice(kgr,3))),"LineWidth",cell2mat(group_choice(kgr,5)),"Color",str2num(cell2mat(group_choice(kgr,3))),'DisplayName',char(group_choice(kgr,2)));
                            else
                                hvol(kgr).plot = plot(xy(:,1),xy(:,2),"LineStyle",char(group_choice(kgr,4)),"Marker",char(group_choice(kgr,6)),"MarkerSize",cell2mat(group_choice(kgr,7)),"LineWidth",cell2mat(group_choice(kgr,5)),"Color",str2num(cell2mat(group_choice(kgr,3))),'DisplayName',char(group_choice(kgr,2)));
                            end
                        else
                            xvalues = [xvalues; xy(:,1)];
                            yvalues = [yvalues; xy(:,2)];
                            tmp = cell(1,length(xy(:,1)));
                            tmp(1,:) = {char(group_choice(kgr,2))};
                            grpvalues = [grpvalues tmp];
                            scol = [scol; str2num(cell2mat(group_choice(kgr,3)))];
                        end

                    end

                end

                if ~plotdist
                    hold off;
                else
                    s=scatterhistogram(xvalues,yvalues,'GroupData',grpvalues,'HistogramDisplayStyle','smooth','LineStyle','-','LegendVisible','on');
                    s.Color  = scol;
                    %s.LegendTitle = 'Group';
                end

            end


            % X-Axis labels
            if row==number_parameter
                par_name = char(parms.Parameter(col));
                unit_name = char(parms.Unit(col));
                if isempty(unit_name)
                    xlabel(par_name)
                else
                    if strcmp(unit_name,'um') || strcmp(unit_name,'micrometer') || strcmp(unit_name,'Micrometer') || strcmp(unit_name,'micrometers') || strcmp(unit_name,'Micrometers')
                        unit_name = '(\mum)';
                    elseif strcmp(unit_name,'um-1') || strcmp(unit_name,'micrometer-1') || strcmp(unit_name,'Micrometer-1') || strcmp(unit_name,'micrometers-1') || strcmp(unit_name,'Micrometers-1')
                        unit_name = '(\mum^{-1})';
                    elseif strcmp(unit_name,'um-2') || strcmp(unit_name,'micrometer-2') || strcmp(unit_name,'Micrometer-2') || strcmp(unit_name,'micrometers-2') || strcmp(unit_name,'Micrometers-2')
                        unit_name = '(\mum^{-2})';   

                    elseif strcmp(unit_name,'nm') || strcmp(unit_name,'nanometer') || strcmp(unit_name,'Nanometer') || strcmp(unit_name,'nanometers') || strcmp(unit_name,'Nanometers')
                        unit_name = '(nm)';
                    elseif strcmp(unit_name,'nm-1') || strcmp(unit_name,'nanometer-1') || strcmp(unit_name,'Nanometer-1') || strcmp(unit_name,'nanometers-1') || strcmp(unit_name,'Nanometers-1')
                        unit_name = '(nm^{-1})';
                    elseif strcmp(unit_name,'nm-2') || strcmp(unit_name,'nanometer-2') || strcmp(unit_name,'Nanometer-2') || strcmp(unit_name,'nanometers-2') || strcmp(unit_name,'Nanometers-2')
                        unit_name = '(nm^{-2})';   

                    else
                        unit_name = ['(' unit_name ')'];
                    end
                    xlabel([par_name ' ' unit_name])
                end
            end
            % Y-Axis labels
            if col==1
                par_name = char(parms.Parameter(row));
                unit_name = char(parms.Unit(row));
                if isempty(unit_name)
                    ylabel(par_name)
                else
                    if strcmp(unit_name,'um') || strcmp(unit_name,'micrometer') || strcmp(unit_name,'Micrometer') || strcmp(unit_name,'micrometers') || strcmp(unit_name,'Micrometers')
                        unit_name = '(\mum)';
                    elseif strcmp(unit_name,'um-1') || strcmp(unit_name,'micrometer-1') || strcmp(unit_name,'Micrometer-1') || strcmp(unit_name,'micrometers-1') || strcmp(unit_name,'Micrometers-1')
                        unit_name = '(\mum^{-1})';
                    elseif strcmp(unit_name,'um-2') || strcmp(unit_name,'micrometer-2') || strcmp(unit_name,'Micrometer-2') || strcmp(unit_name,'micrometers-2') || strcmp(unit_name,'Micrometers-2')
                        unit_name = '(\mum^{-2})';       

                    elseif strcmp(unit_name,'nm') || strcmp(unit_name,'nanometer') || strcmp(unit_name,'Nanometer') || strcmp(unit_name,'nanometers') || strcmp(unit_name,'Nanometers')
                        unit_name = '(nm)';
                    elseif strcmp(unit_name,'nm-1') || strcmp(unit_name,'nanometer-1') || strcmp(unit_name,'Nanometer-1') || strcmp(unit_name,'nanometers-1') || strcmp(unit_name,'Nanometers-1')
                        unit_name = '(nm^{-1})';
                    elseif strcmp(unit_name,'nm-2') || strcmp(unit_name,'nanometer-2') || strcmp(unit_name,'Nanometer-2') || strcmp(unit_name,'nanometers-2') || strcmp(unit_name,'Nanometers-2')
                        unit_name = '(nm^{-2})';  

                    else
                        unit_name = ['(' unit_name ')'];
                    end
                    ylabel([par_name ' ' unit_name])
                end
            end
            set(gca,'Fontname',pformat.fontname,'Fontsize',pformat.size.axes)
            if ~plotdist
                grid(gca,pformat.grid); % Display grid
                set(gca,'XMinorGrid',pformat.minorgrid,'YMinorGrid',pformat.minorgrid);
            end
        end

    end
end

if ~plotdist
    % Legend, shamelessly using eval
    if strcmp(plotper,'volume')
        n = number_volume;
    else
        n = number_group;
    end
    for kvol = 1:1:n
        str_eval = [str_eval 'hvol(' num2str(kvol) ').plot,'];
    end
    str_eval(end) = [];
    str_eval = ['hL = legend([' str_eval ']);'];
    eval(str_eval);
    if strcmp(pformat.legendposition,'South')
        hL.Layout.Tile = 'South';
        hL.NumColumns = pformat.legendcolumn;
    else
        hL.Layout.Tile = 'East';
    end
    hL.FontSize = pformat.size.legend;
end


end