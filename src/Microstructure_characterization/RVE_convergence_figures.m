function [] = RVE_convergence_figures(p)

do_fig_allphase = true;
do_fig_perphase = true;
if p.number_domain_todo==1
    do_fig_allphase = false;
    do_fig_perphase = true;
end

do_fig_allphase = false;

strunit =  p.infovol.unit;
if strcmp(strunit,'um') || strcmp(strunit,'micrometer') || strcmp(strunit,'Micrometer') || strcmp(strunit,'micrometers') || strcmp(strunit,'Micrometers')
    axisunit = '(\mum)';
else
    axisunit = ['(' strunit ')'];
end 

scrsz = get(0,'ScreenSize'); % Screen resolution

p.propertyname(1) = upper(p.propertyname(1)); % Uppercase for first letter

% Title
str_title_1 = ['FOV is cropped along ' p.RVE.RVEconvergence_Crop ': ' p.RVE.RVEconvergence_thatis];
if length(p.Wholevolume_size)==1
    str_title_2 = [char(p.Analysis_name(1)) ' for ' p.propertyname];
else
    str_title_2 = [char(p.Analysis_name(1)) ' and ' char(p.Analysis_name(2)) ' for ' p.propertyname];
end
str_title_3 = p.RVE.choice_RVE;
if ~strcmp(p.RVE.Constantdirection,'n/a')
    str_title_3 = [str_title_3 ': ' p.RVE.Constantdirection];
end
if ~ischar(p.RVE.Aspectratio_RVE)
    str_title_3 = [str_title_3 ', aspect ratio: ' num2str(p.RVE.Aspectratio_RVE)];
end
str_title = {str_title_1,str_title_2,str_title_3};

% xlabel: FOV
str_xlabel={['FOV ' char(p.FOV_length_name(1)) ' ' axisunit]};
if ~isempty( cell2mat(p.FOV_length_name(2)) )
    str_xlabel(2)={['FOV ' char(p.FOV_length_name(2)) ' ' axisunit]};
end

% ylabel: RVE
str = char(p.Analysis_name(1));
idx = find(str == '(');
short = str(idx+1:end-1);
str_ylabel={[short ' ' char(p.RVE_length_name(1)) ' ' axisunit]};
if ~isempty( cell2mat(p.Analysis_name(2)) )
    str = char(p.Analysis_name(2));
    idx = find(str == '(');
    short = str(idx+1:end-1);
    str_ylabel(2)={[short ' ' char(p.RVE_length_name(2)) ' ' axisunit]};
end

str_criterion = '(RelStd_{threshold} = ';
thresholds = p.RVE.threshold_std;
n_threshold = length(thresholds);

current_domain_todo = 0;
for current_domain=1:1:p.number_domain
    if p.todo(current_domain)
        current_domain_todo=current_domain_todo+1;
        domaintodo(current_domain_todo).name = char(p.infovol.phasename(current_domain));
    end
end

%% ONE FIGURE WITH ALL PHASES
if do_fig_allphase

    for kxlabel = 1:length(str_xlabel)
        n_axe = length(str_ylabel);

        % XX1=[]; YY1=[]; ZZ1=[];
        % XX2=[]; YY2=[]; ZZ2=[];
        Fig = figure; % Create figure
        Fig.Color='white'; % Background colour
        Fig.Name= ['Representativity analysis convergence: ' p.propertyname];
        set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
        for id_axe=1:1:n_axe % Iterate over axe
            k_RVEsize = id_axe;
            sub_axes=subplot(1,n_axe,id_axe,'Parent',Fig);
            hold(sub_axes,'on'); % Active subplot
            clear str_legend;
            k_legend=0;
            for current_domain_todo=1:1:p.number_domain_todo
                for k_threshold=1:1:n_threshold
                    % x_ = p.Result_RVEconv(:,1,current_domain_todo,2,k_RVEsize);
                    x_ = p.length_FOV(:,kxlabel);
                    y_ = p.Result_RVEconv(:,k_threshold+1,current_domain_todo,2,k_RVEsize);
                    id0 = find(y_==0);
                    if length(id0)<length(y_)
                        k_legend=k_legend+1;
                        str_legend(k_legend).name = [domaintodo(current_domain_todo).name ' ' str_criterion num2str(thresholds(k_threshold),'%1.1f') '%)'];
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
                        set(h_, 'Color', p.infovol.phasecolor(current_domain_todo,:),'LineWidth',thresholds(k_threshold),'MarkerSize',p.opt.format.axefontsize,'Marker','o');
                    end
                end
            end

            if exist('str_legend','var')
                h_legend = legend(sub_axes,str_legend.name,'Location','best'); % Legend
                h_legend.FontSize = p.opt.format.legendfontsize; % Set fontsize
            end
            % - Axis label
            xlabel(str_xlabel(kxlabel));
            ylabel(str_ylabel(k_RVEsize));
            grid(sub_axes,p.opt.format.grid); % Display grid
            set(sub_axes,'XMinorGrid',p.opt.format.minorgrid,'YMinorGrid',p.opt.format.minorgrid); % Display grid for minor thicks
            set(sub_axes,'FontName',p.opt.format.fontname,'FontSize',p.opt.format.axefontsize); % Fontname and fontsize
            h_title.FontSize = p.opt.format.titlefontsize; % Set title fontsize
            % axis equal;
            hold(sub_axes,'off'); % Relase figure
        end
        sgtitle(Fig,str_title,'FontWeight','bold','FontSize',p.opt.format.sgtitlefontsize,'FontName',p.opt.format.fontname);
        if p.opt.save.savefig % Save figure
            filename= [p.propertyname '_' p.RVE.savename '_convergenceAll_' num2str(kxlabel)];
            function_savefig(Fig, p.savefolder, filename, p.opt.save); % Call function
        end
        if p.opt.format.autoclosefig
            close(Fig); % Do not keep open figures
        end
    end
end

%% ONE FIGURE PER PHASES
if do_fig_perphase

    for kxlabel = 1:length(str_xlabel)
        r_axe = length(str_ylabel);
        c_axe = p.number_domain_todo;

        Fig = figure; % Create figure
        Fig.Color='white'; % Background colour
        Fig.Name= ['Representativity analysis convergence: ' p.propertyname];

        val = p.Result_RVEconv(:,2:end,:,2,:);
        idx = find(val>0);
        if ~isempty(idx)
            [~, I2, ~, ~, ~] = ind2sub(size(val),idx);
            ploted_threshold = thresholds(unique(I2));
            col=turbo(256);
            vq_col = round(interp1([min(ploted_threshold) max(ploted_threshold)],[256 1],ploted_threshold,'linear'));
        end

        Fig.Color='white'; % Background colour
        set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
        id_axe=0;
        for r=1:1:r_axe
            k_RVEsize = r;
            for c=1:1:c_axe
                current_domain_todo=c;
                id_axe=id_axe+1;
                sub_axes=subplot(r_axe,c_axe,id_axe,'Parent',Fig);
                hold(sub_axes,'on'); % Active subplot

                clear str_legend;
                k_legend=0;
                for k_threshold=1:1:n_threshold
                    %x_ = p.Result_RVEconv(:,1,current_domain_todo,2,k_conv);
                    x_ = p.length_FOV(:,kxlabel);
                    y_ = p.Result_RVEconv(:,k_threshold+1,current_domain_todo,2,k_RVEsize);

                    id0 = find(y_==0);
                    if length(id0)<length(y_)
                        k_legend=k_legend+1;
                        str_legend(k_legend).name = [domaintodo(current_domain_todo).name ' ' str_criterion num2str(thresholds(k_threshold),'%1.1f') '%)'];
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
                xlabel(str_xlabel(kxlabel));
                ylabel(str_ylabel(k_RVEsize));
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
            filename= [p.propertyname '_' p.RVE.savename '_convergence_' num2str(kxlabel)];
            function_savefig(Fig, p.savefolder, filename, p.opt.save); % Call function
        end
        if p.opt.format.autoclosefig
            close(Fig); % Do not keep open figures
        end
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
