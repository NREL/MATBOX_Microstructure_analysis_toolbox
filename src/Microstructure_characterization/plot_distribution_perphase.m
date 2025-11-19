function [ ] = plot_distribution_perphase(LabelsMap,phaselabel,phasename,phasecolor,Array,stats,p,optformat,optsave,folder,filename)

%% STATISTICS
nlabel = length(phaselabel);
if isempty(stats)
    stats = zeros(nlabel,6);
    for k=1:nlabel
        BWlog = LabelsMap==phaselabel(k);
        vals = double(Array(BWlog));
        vals = reshape(vals,[1 numel(vals)]);
        if numunique(vals)==1 % Std can provide non-zero (numerical error) if vals is a long array with the same repeating value
            stats(k,:) = [vals(1) vals(1) vals(1) 0 vals(1) vals(1)];
        else
            stats(k,:) = [mean(vals) median(vals) mode(vals) std(vals) min(vals) max(vals)];
        end
    end
end

%% PDF/CDF
if max(stats(:,4))==0 % Standard deviation is zero
    p.plots = {'bar'}; % Overwritte
end
if ismember('pdf',p.plots) || ismember('cdf',p.plots)
    pars.round_value = p.round_value;
    pars.smooth_cumulative_fct = p.smooth_cumulative_fct;
   
    sheet=0;
    for k=1:nlabel
        idx = find(LabelsMap==phaselabel(k));
        vals = Array(idx);
        if length(unique(vals))>1
            [pd(k).res, ~] = Function_probability_density(vals,[],pars);
            
            sheet=sheet+1;
            varnames = [{[p.array_name ' ' p.array_unit]} {'cdf'}];
            DATA_writetable.sheet(sheet).name = [char(phasename(k)) ' cdf'];
            DATA_writetable.sheet(sheet).table = array2table([pd(k).res.cumulative_fct(:,1) pd(k).res.cumulative_fct(:,2)],"VariableNames",varnames);
            if ~isempty(pd(k).res.smoothed_cumulative_fct)
                sheet=sheet+1;
                varnames = [{[p.array_name ' ' p.array_unit]} {'cdf (smoothed)'}];
                DATA_writetable.sheet(sheet).name = [char(phasename(k)) ' cdf (smoothed)'];
                DATA_writetable.sheet(sheet).table = array2table([pd(k).res.smoothed_cumulative_fct(:,1) pd(k).res.smoothed_cumulative_fct(:,2)],"VariableNames",varnames);
            end
            sheet=sheet+1;
            varnames = [{[p.array_name ' ' p.array_unit]} {'pdf'}];
            DATA_writetable.sheet(sheet).name = [char(phasename(k)) ' pdf'];
            DATA_writetable.sheet(sheet).table = array2table([pd(k).res.probability_density_fct(:,1) pd(k).res.probability_density_fct(:,2)],"VariableNames",varnames);
            if ~isempty(pd(k).res.smoothed_probability_density_fct)
                sheet=sheet+1;
                varnames = [{[p.array_name ' ' p.array_unit]} {'pdf (smoothed)'}];
                DATA_writetable.sheet(sheet).name = [char(phasename(k)) ' pdf (smoothed)'];
                DATA_writetable.sheet(sheet).table = array2table([pd(k).res.smoothed_probability_density_fct(:,1) pd(k).res.smoothed_probability_density_fct(:,2)],"VariableNames",varnames);
            end
        else
            pd(k).res=[];
        end
    end
    if sheet>0
        Function_Writetable(folder,[filename '_densityfct'],DATA_writetable)
    end
end

%% FIGURE
fig = figure;
fig.Color = 'w';
fig.Name = p.array_name;

nplot = length(p.plots);
tiledlayout(1,nplot,'TileSpacing',optformat.tile_spacing,'Padding',optformat.layout_padding);

h1 = optsave.Height_foroneplot;
w1 = h1;
optsave.Width = nplot*w1; % + (nplot-1) * WidthBtw;
optsave.Height = h1;

for k=1:nplot
    ax = nexttile;
    if strcmp(char(p.plots(k)),'bar')
        b=bar(ax,phasename,stats(:,1),1.0);
        b.FaceColor = 'flat';
        b.Labels = round(stats(:,1),3);
        for current_phase=1:1:nlabel
            b.CData(current_phase,:)=phasecolor(current_phase,:);
        end  
        xlabel(ax,"Phase")
        ylabel(ax,[p.array_name ', average per phase'])
        ysecondarylabel(ax,p.array_unit)
        ax.Box = "off";

    elseif strcmp(char(p.plots(k)),'bar(mean and std)')
        tmpname = cell(2*nlabel,1);
        tmpstat = zeros(2*nlabel,1);
        tmpcolor = zeros(2*nlabel,3);
        for kk=1:nlabel
            tmpname(2*kk-1) = {[char(phasename(kk)) '(mean)']};
            tmpname(2*kk) = {[char(phasename(kk)) '(std)']};
            tmpstat(2*kk-1) = stats(kk,1);
            tmpstat(2*kk) = stats(kk,4);
            tmpcolor(2*kk-1,:) = phasecolor(kk,:);
            tmpcolor(2*kk,:) = (1+phasecolor(kk,:))/2;
        end
        b=bar(ax,tmpname,tmpstat,1.0);
        b.FaceColor = 'flat';
        b.Labels = round(tmpstat,3);
        for kk=1:2*nlabel
            b.CData(kk,:)=tmpcolor(kk,:);
        end  
        xlabel(ax,"Phase")
        ylabel(ax,p.array_name)
        ysecondarylabel(ax,p.array_unit)
        ax.Box = "off";

        % b=bar(ax,phasename,stats(:,1:2),1.0);
        % b(1).Labels = round(stats(:,1),3);
        % b(2).Labels = round(stats(:,2),3);
        % set(b(1),'DisplayName','Mean','FaceColor',[0 0 0]);
        % set(b(2),'DisplayName','Std','FaceColor',[0.5 0.5 0.5]);
        % hl = legend;
        % hl.FontSize = optformat.legendfontsize;
        % hl.Location = 'best';
        % hl.IconColumnWidth = 20;
        % hl.LineWidth = 1;
        % xlabel(ax,"Phase")
        % ylabel(ax,p.array_name)  
        % ysecondarylabel(ax,p.array_unit)
        % ax.Box = "off";

    elseif strcmp(char(p.plots(k)),'pdf')
        hold(ax,'on');
        for k=1:nlabel
            if ~isempty(pd(k).res)
                if pars.smooth_cumulative_fct
                    plot(pd(k).res.probability_density_fct(:,1), pd(k).res.probability_density_fct(:,2),'Color',(1+phasecolor(k,:))/2,'DisplayName',[char(phasename(k)) ', x_{50} (median) = ' num2str(pd(k).res.x50) ' ' p.array_unit],'LineWidth',optformat.linewidth/2);
                    plot(pd(k).res.smoothed_probability_density_fct(:,1), pd(k).res.smoothed_probability_density_fct(:,2),'Color',phasecolor(k,:),'DisplayName',[char(phasename(k)) ' (smoothed), x_{50} (median) = ' num2str(pd(k).res.smoothed_x50) ' ' p.array_unit],'LineWidth',optformat.linewidth,'LineStyle','-');
                else
                    plot(pd(k).res.probability_density_fct(:,1), pd(k).res.probability_density_fct(:,2),'Color',phasecolor(k,:),'DisplayName',[char(phasename(k)) ', x_{50} (median) = ' num2str(pd(k).res.x50) ' ' p.array_unit],'LineWidth',optformat.linewidth);
                end              
            end
        end
        hold(ax,'off');
        xlabel(ax,p.array_name)
        xsecondarylabel(ax,p.array_unit)
        ylabel(ax,"Probability density function (a.u)")
        hl = legend;
        hl.FontSize = optformat.legendfontsize;
        hl.Location = 'best';
        hl.IconColumnWidth = 20;
        hl.LineWidth = 1;


    elseif strcmp(char(p.plots(k)),'cdf')
        hold(ax,'on');
        for k=1:nlabel
            if ~isempty(pd(k).res)
                if pars.smooth_cumulative_fct
                    plot(pd(k).res.cumulative_fct(:,1), pd(k).res.cumulative_fct(:,2),'Color',(1+phasecolor(k,:))/2,'DisplayName',[char(phasename(k)) ', x_{50} (median) = ' num2str(pd(k).res.x50) ' ' p.array_unit],'LineWidth',optformat.linewidth/2);
                    plot(pd(k).res.smoothed_cumulative_fct(:,1), pd(k).res.smoothed_cumulative_fct(:,2),'Color',phasecolor(k,:),'DisplayName',[char(phasename(k)) ' (smoothed), x_{50} (median) = ' num2str(pd(k).res.smoothed_x50) ' ' p.array_unit],'LineWidth',optformat.linewidth,'LineStyle','-');
                else
                    plot(pd(k).res.cumulative_fct(:,1), pd(k).res.cumulative_fct(:,2),'Color',phasecolor(k,:),'DisplayName',[char(phasename(k)) ', x_{50} (median) = ' num2str(pd(k).res.x50) ' ' p.array_unit],'LineWidth',optformat.linewidth);
                end
            end
        end  
        n=0;
        for k=1:nlabel
            if ~isempty(pd(k).res)
                if pars.smooth_cumulative_fct
                    plot([0 pd(k).res.x50 pd(k).res.x50],[0.5 0.5 0.0],'Color',(1+phasecolor(k,:))/2,'LineWidth',1,'LineStyle','--');
                    plot([0 pd(k).res.smoothed_x50 pd(k).res.smoothed_x50],[0.5 0.5 0.0],'Color',phasecolor(k,:),'LineWidth',1,'LineStyle','--');
                    n=n+2;
                else
                    plot([0 pd(k).res.x50 pd(k).res.x50],[0.5 0.5 0.0],'Color',phasecolor(k,:),'LineWidth',1,'LineStyle','--');
                    n=n+1;
                end
            end
        end  
        ylim(ax,[0,1]);
        hold(ax,'off');
        xlabel(ax,p.array_name)
        xsecondarylabel(ax,p.array_unit)
        ylabel(ax,"Cumulative density function")
        hl = legend;
        hl.String(end-n+1:end)=[];
        hl.FontSize = optformat.legendfontsize;
        hl.Location = 'best';
        hl.IconColumnWidth = 20;
        hl.LineWidth = 1;
    end

    % Common
    ax.XLabel.FontWeight = "bold";
    ax.YLabel.FontWeight = "bold";
    ax.LineWidth = optformat.linewidth;
    ax.FontName = optformat.fontname;
    ax.FontSize = optformat.axefontsize;
    if optformat.grid
        grid(ax,'on');
        set(ax,'XMinorGrid',optformat.minorgrid,'YMinorGrid',optformat.minorgrid,'GridLineWidth',1);
    end

end

if optformat.includefilenameintitle
    if nplot==1
        title(p.inputfilename,'Interpreter','none','FontSize',optformat.titlefontsize);
    else
        sgtitle(p.inputfilename,'Interpreter','none','FontSize',optformat.sgtitlefontsize);
    end
end

%% SAVE FIGURE
if optsave.savefig % Save figure
    function_savefig(fig, folder, filename, optsave); % Call function
end
if optsave.autoclosefig
    close(fig); % Do not keep open figures
end

end