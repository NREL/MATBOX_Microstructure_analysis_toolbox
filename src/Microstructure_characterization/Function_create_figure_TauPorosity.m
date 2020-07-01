function [] = Function_create_figure_TauPorosity(p,datafrom)

if strcmp(datafrom, 'Voxel_size_dependence')

    Fig = figure; % Create figure
    Fig.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)/2]); % Full screen figure    
    % - Create axes as a subplot
    for id_axe=1:1:3 % id_axe = current_direction
        clear str_legend
        sub_axes=subplot(1,3,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        h_title=title (p.INFO.direction(id_axe).name); % Title
        k_legend=0;
        for current_phase=1:1:p.number_phase % Loop over phases
            % Select data
            tau = p.direction(id_axe).property_voxelsizedependence(2:end,current_phase+1,1);
            epsilon = p.direction(id_axe).property_voxelsizedependence(2:end,current_phase+1,4);
            % Calculate tortuosity porosity coefficients
            [results] = Function_TortuosityPorosity_fitting(epsilon,tau);
            if results.success
                h_point=plot(epsilon,tau);
                set(h_point, 'Color', p.INFO.phase(current_phase).color,'MarkerSize',p.OPTIONS.Fontsize_axe,'Marker','o','LineWidth',p.OPTIONS.Linewidth,'Linestyle','none');
                k_legend=k_legend+1;
                str_legend(k_legend).name = [p.INFO.phase(current_phase).name ', points'];
                if isfield(results,'epsilon_fit') && isfield(results,'tau_fit') && isfield(results,'gamma') && isfield(results,'alpha')
                    h_fit=plot(results.epsilon_fit,results.tau_fit);
                    set(h_fit, 'Color', p.INFO.phase(current_phase).color,'LineWidth',p.OPTIONS.Linewidth,'Linestyle','--');
                    k_legend=k_legend+1;
                    if results.gamma>0
                        str_legend(k_legend).name = [p.INFO.phase(current_phase).name ', fit: \tau = ' num2str(results.alpha,'%1.3f') ' \times \epsilon power(-' num2str(results.gamma,'%1.3f') ')'];
                    else
                        str_legend(k_legend).name = [p.INFO.phase(current_phase).name ', fit: \tau = ' num2str(results.alpha,'%1.3f') ' \times \epsilon power(' num2str(-results.gamma,'%1.3f') ')'];
                    end
                end
            end
        end
        % Axis label
        xlabel('Volume fraction (\epsilon)');
        ylabel('Tortuosity factor (\tau)');
        % Legend
        h_legend = legend(sub_axes,str_legend.name,'Location','best');
        % - Grid
        if strcmp(p.OPTIONS.grid,'on')
            grid(sub_axes,'on'); % Display grid
            set(sub_axes,'XMinorGrid',p.OPTIONS.minorgrid,'YMinorGrid',p.OPTIONS.minorgrid); % Display grid for minor thicks
        end
        set(sub_axes,'FontName',p.OPTIONS.fontname,'FontSize',p.OPTIONS.Fontsize_axe); % Fontname and fontsize
        h_title.FontSize = p.OPTIONS.Fontsize_title; % Set title fontsize
        h_legend.FontSize = p.OPTIONS.Fontsize_legend; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
    sgtitle(Fig,'Tortuosity factor - volume fraction correlation, from voxel size dependence analysis','FontWeight','bold','FontSize',p.OPTIONS.Fontsize_title+2,'FontName',p.OPTIONS.fontname);

else
    
end

if p.OPTIONS.save_fig == true % Save figure
    function_savefig(Fig, p.Current_folder, p.filename, p.OPTIONS); % Call function
end
if p.OPTIONS.closefigureaftercreation == true
    close(Fig); % Do not keep open figures
end

end

