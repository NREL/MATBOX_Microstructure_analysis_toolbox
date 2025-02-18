function [] = function_fractalfig(propertyname,filename, Current_folder, box_lengths,N,fractal_dimension_convergence,Topology_dimension,number_phase,p,opt,infovol)

scrsz = get(0,'ScreenSize'); % Screen resolution
propertyname(1) = upper(propertyname(1));
Fig = figure; % Create figure
Fig.Name= [propertyname ' fractal dimension (box counting)']; % Figure name
Fig.Color='white'; % Background colour
set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*4/5 scrsz(4)*4/5]); % Full screen figure

for id_axe=1:1:4
    sub_axes=subplot(2,2,id_axe,'Parent',Fig);
    hold(sub_axes,'on'); % Active subplot
    if id_axe==1
        sub_axes=subplot(2,2,id_axe,'Parent',Fig,'XScale', 'log', 'YScale', 'log');
        h_title=title ('Fractal dimension fit (1/2)'); % Set title font

        x=1./box_lengths;
        current_phase_todo = 0;
        for current_phase=1:1:number_phase % Loop over phases
            if p.todo(current_phase)
                current_phase_todo=current_phase_todo+1;
                y=N(:,current_phase_todo);
                plot(x,y,'Color', infovol.phasecolor(current_phase,:),'LineWidth',opt.format.linewidth,'DisplayName',char(infovol.phasename(current_phase,1)));
            end
        end
        xlabel('1 / box length (logarithmic scale)');
        ylabel('N (logarithmic scale)');
        % axis 'equal';

    elseif id_axe==2
        sub_axes=subplot(2,2,id_axe,'Parent',Fig);
        h_title=title ('Fractal dimension fit (2/2)'); % Set title font
        current_phase_todo = 0;
        for current_phase=1:1:number_phase % Loop over phases
            if p.todo(current_phase)
                current_phase_todo=current_phase_todo+1;
                plot(fractal_dimension_convergence(:,1,current_phase_todo),fractal_dimension_convergence(:,3,current_phase_todo),'Color', infovol.phasecolor(current_phase,:),'LineWidth',opt.format.linewidth,'DisplayName',char(infovol.phasename(current_phase,1)));
            end
        end
        xlabel('Number of box length used to fit the fractal dimension (from small to large length)');
        ylabel('Polyfit residual norm');

    elseif id_axe==3
        sub_axes=subplot(2,2,id_axe,'Parent',Fig);
        h_title=title ('Fractal dimension'); % Set title font
        current_phase_todo = 0;
        for current_phase=1:1:number_phase % Loop over phases
            if p.todo(current_phase)
                current_phase_todo=current_phase_todo+1;
                plot(fractal_dimension_convergence(:,1,current_phase_todo),fractal_dimension_convergence(:,2,current_phase_todo),'Color', infovol.phasecolor(current_phase,:),'LineWidth',opt.format.linewidth,'DisplayName',char(infovol.phasename(current_phase,1)));
            end
        end
        plot([min(fractal_dimension_convergence(:,1,1)) max(fractal_dimension_convergence(:,1,1))],[Topology_dimension Topology_dimension],'Color','k','LineWidth',1,'LineStyle','--','DisplayName','Topology dimension');
        xlabel('Number of box length used to fit the fractal dimension (from small to large length)');
        ylabel('Fractal dimension');

    elseif id_axe==4
        sub_axes=subplot(2,2,id_axe,'Parent',Fig);
        h_title=title ('Fractal propensity'); % Set title font
        current_phase_todo = 0;
        for current_phase=1:1:number_phase % Loop over phases
            if p.todo(current_phase)
                current_phase_todo=current_phase_todo+1;
                plot(fractal_dimension_convergence(:,1,current_phase_todo),abs(Topology_dimension-fractal_dimension_convergence(:,2,current_phase_todo)),'Color', infovol.phasecolor(current_phase,:),'LineWidth',opt.format.linewidth,'DisplayName',char(infovol.phasename(current_phase,1)));
            end
        end
        plot([min(fractal_dimension_convergence(:,1,1)) max(fractal_dimension_convergence(:,1,1))],[0 0],'Color','k','LineWidth',1,'LineStyle','--','DisplayName','Zero-fractal');
        xlabel('Number of box length used to fit the fractal dimension (from small to large length)');
        ylabel('abs(Topology dimension - fractal dimension)');
    end
    % Legend
    h_legend = legend(sub_axes,'Location','best');
    % - Grid
    grid(sub_axes,'on'); % Display grid
    set(sub_axes,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks
    set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
    h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
    h_legend.FontSize = opt.format.legendfontsize; % Set title fontsize
    hold(sub_axes,'off'); % Relase figure
end
sgtitle(Fig,Fig.Name,'FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
if opt.save.savefig % Save figure
    function_savefig(Fig, Current_folder, filename, opt.save); % Call function
end
if opt.format.autoclosefig
    close(Fig); % Do not keep open figures
end

end