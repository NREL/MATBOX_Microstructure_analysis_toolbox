function [] = Function_create_figures_nestedRVE2(p)

do_fig_allphase = true;
do_fig_perphase = true;
if p.number_phase_todo==1
    do_fig_allphase = false;
    do_fig_perphase = true;
end

strunit =  p.infovol.unit;
if strcmp(strunit,'um') || strcmp(strunit,'micrometer') || strcmp(strunit,'Micrometer') || strcmp(strunit,'micrometers') || strcmp(strunit,'Micrometers')
    axisunit = '(\mum)';
else
    axisunit = ['(' strunit ')'];
end 

scrsz = get(0,'ScreenSize'); % Screen resolution

p.propertyname(1) = upper(p.propertyname(1)); % Uppercase for first letter

Fig = figure; % Create figure
if strcmp(p.RVE.type,'A') || strcmp(p.RVE.type,'B') || strcmp(p.RVE.type,'C') || strcmp(p.RVE.type,'D')
    Fig.Name= ['Representative volume element convergence analysis: ' p.propertyname];
else
    Fig.Name= ['Convergence analysis:' p.propertyname];
end

if strcmp(p.RVE.type,'A') || strcmp(p.RVE.type,'B')
    str_title = {Fig.Name,p.RVE.choice_RVE,['Aspect ratio: ' num2str(p.RVE.Aspectratio_RVE)]};
elseif strcmp(p.RVE.type,'C') || strcmp(p.RVE.type,'D') 
    str_title = {Fig.Name,p.RVE.choice_RVE};
elseif strcmp(p.RVE.type,'E') || strcmp(p.RVE.type,'F')
    str_title = {Fig.Name,p.RVE.choice_Onsub,['Aspect ratio: ' num2str(p.RVE.Aspectratio_Onesub)]};
elseif strcmp(p.RVE.type,'G')
    str_title = {Fig.Name,p.RVE.choice_Onsub,p.RVE.Growthdirection};
elseif strcmp(p.RVE.type,'H')
    str_title = {Fig.Name,p.RVE.choice_Onsub,p.RVE.Constantdirection};    
end

if p.dimension==3
    str_xlabel={['FOV cubic root ' axisunit], []};
    if strcmp(p.RVE.type,'C') || strcmp(p.RVE.type,'H')
        str_xlabel(2)={['FOV square root ' axisunit]};
        ky = 2;
    elseif strcmp(p.RVE.type,'D')
        str_xlabel(2)={['FOV length ' axisunit]};
        ky = 3;
    elseif strcmp(p.RVE.type,'F')
        str_xlabel(2)={['FOV length ' axisunit]};
        ky = 3;
    else
        foo=1;
    end

else
    str_xlabel={['FOV square root ' axisunit], []};
    if strcmp(p.RVE.type,'C') || strcmp(p.RVE.type,'H')
        str_xlabel(2)={['FOV length ' axisunit]};
        ky = 2;
    elseif strcmp(p.RVE.type,'D')
        str_xlabel(2)={['FOV length ' axisunit]};
        ky = 3;
    elseif strcmp(p.RVE.type,'F')
        str_xlabel(2)={['FOV length ' axisunit]};
        ky = 3;
    else
        foo=1;
    end
end

if p.dimension == 3
    if strcmp(p.RVE.type,'A') || strcmp(p.RVE.type,'B')
        str_ylabel={['RVE cubic root ' axisunit], []};
    elseif strcmp(p.RVE.type,'C')
        str_ylabel={['RVE cubic root ' axisunit], ['RSA square root ' axisunit]};
    elseif strcmp(p.RVE.type,'D')
        str_ylabel={['RVE cubic root ' axisunit], ['Representative length ' axisunit]};
    elseif strcmp(p.RVE.type,'E') || strcmp(p.RVE.type,'F')
        str_ylabel={['Convergence volume, cubic root ' axisunit], []};
    elseif strcmp(p.RVE.type,'G')
        str_ylabel={['Convergence volume, cubic root ' axisunit], ['Convergence volume, length ' axisunit]};
    elseif strcmp(p.RVE.type,'H')
        str_ylabel={['Convergence volume, cubic root ' axisunit], ['Convergence volume, square root ' axisunit]};
    end
else
    if strcmp(p.RVE.type,'A') || strcmp(p.RVE.type,'B')
        str_ylabel={['RSA square root ' axisunit], []};
    elseif strcmp(p.RVE.type,'C')
        str_ylabel={['RSA square root ' axisunit], ['Representative length ' axisunit]};
    elseif strcmp(p.RVE.type,'D')
        str_ylabel={['RSA square root ' axisunit], ['Representative length ' axisunit]};
    elseif strcmp(p.RVE.type,'E') || strcmp(p.RVE.type,'F')
        str_ylabel={['Convergence area, square root ' axisunit], []};
    elseif strcmp(p.RVE.type,'G')
        str_ylabel={['Convergence area, square root ' axisunit], ['Convergence area, length ' axisunit]};
    elseif strcmp(p.RVE.type,'H')
        str_ylabel={['Convergence area, square root ' axisunit], ['Convergence area, length ' axisunit]};
    end


end


if strcmp(p.RVE.type,'A') || strcmp(p.RVE.type,'B') || strcmp(p.RVE.type,'C') || strcmp(p.RVE.type,'D')
    str_criterion = '(RelStd_{threshold} = ';
    thresholds = p.RVE.threshold_std;
else
    thresholds = p.RVE.threshold_reldiff;
    str_criterion = '(RelDiff_{threshold} = ';
end
n_threshold = length(thresholds);

current_phase_todo = 0;
for current_phase=1:1:p.number_phase
    if p.todo(current_phase)
        current_phase_todo=current_phase_todo+1;
        phasetodo(current_phase_todo).name = char(p.infovol.phasename(current_phase));
    end
end

%% ONE FIGURE WITH ALL PHASES
if do_fig_allphase
    if strcmp(p.RVE.type,'C') || strcmp(p.RVE.type,'H')
        n_axe = 2;
        ky = 2;
    elseif strcmp(p.RVE.type,'D')
        n_axe = 2;
        ky = 3;        
    elseif strcmp(p.RVE.type,'G')
        n_axe = 2;
        ky = 3;
    else
        n_axe = 1;
    end

    % XX1=[]; YY1=[]; ZZ1=[];
    % XX2=[]; YY2=[]; ZZ2=[];
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
    for id_axe=1:1:n_axe % Iterate over axe
        sub_axes=subplot(1,n_axe,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        clear str_legend;
        k_legend=0;
        for current_phase_todo=1:1:p.number_phase_todo
            for k_threshold=1:1:n_threshold

                if id_axe==1
                    x_ = p.Result_nestedRVE(:,1,current_phase_todo,2,1);
                    y_ = p.Result_nestedRVE(:,k_threshold+1,current_phase_todo,2,1);
                else
                    x_ = p.Result_nestedRVE(:,1,current_phase_todo,2,ky);
                    y_ = p.Result_nestedRVE(:,k_threshold+1,current_phase_todo,2,ky);
                end
                id0 = find(y_==0);
                if length(id0)<length(y_)
                    k_legend=k_legend+1;
                    str_legend(k_legend).name = [phasetodo(current_phase_todo).name ' ' str_criterion num2str(thresholds(k_threshold),'%1.1f') '%)'];
                    x_(id0)=NaN; % Do not display
                    y_(id0)=NaN;
                    h_=plot(x_,y_);
                    % x_(id0)=[];
                    % y_(id0)=[];
                    %if id_axe==1
                    %    XX1 =[XX1; x_]; YY1 =[YY1; y_]; ZZ1 =[ZZ1; ones(length(x_),1)*thresholds(k_threshold)];
                    %else
                    %    XX2 =[XX2; x_]; YY2 =[YY2; y_]; ZZ2 =[ZZ2; ones(length(x_),1)*thresholds(k_threshold)];
                    %end
                    % Colors, thickness, markers
                    set(h_, 'Color', p.infovol.phasecolor(current_phase_todo,:),'LineWidth',thresholds(k_threshold),'MarkerSize',p.opt.format.axefontsize,'Marker','o');
                end
            end
        end

        if exist('str_legend','var')
            h_legend = legend(sub_axes,str_legend.name,'Location','best'); % Legend
            h_legend.FontSize = p.opt.format.legendfontsize; % Set fontsize
        end
        % - Axis label
        xlabel(str_xlabel(id_axe));
        ylabel(str_ylabel(id_axe));
        grid(sub_axes,p.opt.format.grid); % Display grid
        set(sub_axes,'XMinorGrid',p.opt.format.minorgrid,'YMinorGrid',p.opt.format.minorgrid); % Display grid for minor thicks
        set(sub_axes,'FontName',p.opt.format.fontname,'FontSize',p.opt.format.axefontsize); % Fontname and fontsize
        h_title.FontSize = p.opt.format.titlefontsize; % Set title fontsize
        % axis equal;
        hold(sub_axes,'off'); % Relase figure
    end
    sgtitle(Fig,str_title,'FontWeight','bold','FontSize',p.opt.format.sgtitlefontsize,'FontName',p.opt.format.fontname);
    if p.opt.save.savefig % Save figure
        filename= [p.propertyname '_' p.RVE.savename '_nested'];
        function_savefig(Fig, p.savefolder, filename, p.opt.save); % Call function
    end
    if p.opt.format.autoclosefig
        close(Fig); % Do not keep open figures
    end
end

%% ONE FIGURE PER PHASES
if do_fig_perphase

    Fig = figure; % Create figure
    if strcmp(p.RVE.type,'A') || strcmp(p.RVE.type,'B') || strcmp(p.RVE.type,'C') || strcmp(p.RVE.type,'D')
        Fig.Name= ['Representative volume element convergence analysis: ' p.propertyname];
    else
        Fig.Name= ['Convergence analysis:' p.propertyname];
    end

    val = p.Result_nestedRVE(:,2:end,:,2,:);
    idx = find(val>0);
    if ~isempty(idx)
        [~, I2, ~, ~, ~] = ind2sub(size(val),idx);
        ploted_threshold = thresholds(unique(I2));
        col=turbo(256);
        vq_col = round(interp1([min(ploted_threshold) max(ploted_threshold)],[256 1],ploted_threshold,'linear'));
    end

    if strcmp(p.RVE.type,'C') || strcmp(p.RVE.type,'H')
        r_axe = 2;
        ky = 2;
    elseif strcmp(p.RVE.type,'D')
        r_axe = 2;
        ky = 3;        
    elseif strcmp(p.RVE.type,'G')
        r_axe = 2;
        ky = 3;
    else
        r_axe = 1;
    end
    c_axe = p.number_phase_todo;
    n_axe = r_axe*c_axe;

    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
    id_axe=0;
    for r=1:1:r_axe
        for c=1:1:c_axe
            current_phase_todo=c;
            id_axe=id_axe+1;
            sub_axes=subplot(r_axe,c_axe,id_axe,'Parent',Fig);
            hold(sub_axes,'on'); % Active subplot

            clear str_legend;
            k_legend=0;
            for k_threshold=1:1:n_threshold
                if r==1
                    x_ = p.Result_nestedRVE(:,1,current_phase_todo,2,1);
                    y_ = p.Result_nestedRVE(:,k_threshold+1,current_phase_todo,2,1);
                else
                    x_ = p.Result_nestedRVE(:,1,current_phase_todo,2,ky);
                    y_ = p.Result_nestedRVE(:,k_threshold+1,current_phase_todo,2,ky);
                end
                id0 = find(y_==0);
                if length(id0)<length(y_)
                    k_legend=k_legend+1;
                    str_legend(k_legend).name = [phasetodo(current_phase_todo).name ' ' str_criterion num2str(thresholds(k_threshold),'%1.1f') '%)'];
                    x_(id0)=NaN; % Do not display
                    y_(id0)=NaN;
                    %x_(id0)=[]; % Do not display
                    %y_(id0)=[];                    
                    h_=plot(x_,y_);
                    % Colors, thickness, markers
                    idc = find(ploted_threshold==thresholds(k_threshold));
                    set(h_, 'Color',col(vq_col(idc),:), 'LineWidth',p.opt.format.linewidth,'MarkerSize',p.opt.format.axefontsize,'Marker','o');
                end
            end
            if exist('str_legend','var')
                h_legend = legend(sub_axes,str_legend.name,'Location','best'); % Legend
                h_legend.FontSize = p.opt.format.legendfontsize; % Set fontsize
            end
            % - Axis label
            xlabel(str_xlabel(r));
            ylabel(str_ylabel(r));
            grid(sub_axes,p.opt.format.grid); % Display grid
            set(sub_axes,'XMinorGrid',p.opt.format.minorgrid,'YMinorGrid',p.opt.format.minorgrid); % Display grid for minor thicks
            set(sub_axes,'FontName',p.opt.format.fontname,'FontSize',p.opt.format.axefontsize); % Fontname and fontsize
            h_title.FontSize = p.opt.format.titlefontsize; % Set title fontsize
            % axis equal;
            hold(sub_axes,'off'); % Relase figure
        end
    end
        
    sgtitle(Fig,str_title,'FontWeight','bold','FontSize',p.opt.format.sgtitlefontsize,'FontName',p.opt.format.fontname);
    if p.opt.save.savefig % Save figure
        filename= [p.propertyname '_' p.RVE.savename '_nested2'];
        function_savefig(Fig, p.savefolder, filename, p.opt.save); % Call function
    end
    if p.opt.format.autoclosefig
        close(Fig); % Do not keep open figures
    end
end

%% CONTOUR MAP

% F = scatteredInterpolant(XX1,YY1,ZZ1,'linear','none');
% nx=100; ny=100;
% x = linspace(min(XX1),max(XX2),nx);
% y = linspace(min(YY1),max(YY2),ny);
% [Xq,Yq] = meshgrid(x,y);
% Vq = F(Xq,Yq);
% Fig=figure
% %pcolor(Xq,Yq,Vq) 
% contourf(Xq,Yq,Vq,10) 

end
