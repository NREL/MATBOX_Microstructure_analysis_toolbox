function [Fig_] = function_2Dmap_3Darray(Microstructure)

domain_size=size(Microstructure);
dimension=length(domain_size);
if dimension==3
    
    %% CALCULATION
    
    % Normal to axe 1
    greylevel_normal(1).map =  zeros(domain_size(2),domain_size(3));
    for x=1:1:domain_size(2)
        for y=1:1:domain_size(3)
            line=reshape(Microstructure(:,x,y),[domain_size(1), 1]);
            greylevel_normal(1).map(x,y)=mean(line);
        end
    end
    
    % Normal to axe 2
    greylevel_normal(2).map =  zeros(domain_size(1),domain_size(3));
    for x=1:1:domain_size(1)
        for y=1:1:domain_size(3)
            line=reshape(Microstructure(x,:,y),[domain_size(2), 1]);
            greylevel_normal(2).map(x,y)=mean(line);
        end
    end
    
    % Normal to axe 3
    greylevel_normal(3).map =  zeros(domain_size(1),domain_size(2));
    for x=1:1:domain_size(1)
        for y=1:1:domain_size(2)
            line=reshape(Microstructure(x,y,:),[domain_size(3), 1]);
            greylevel_normal(3).map(x,y)=mean(line);
        end
    end
    
    %% FIGURE
    Fig_ = figure;
    Fig_.Name= 'Mean gray level value';
    Fig_.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',[scrsz(1) scrsz(2) scrsz(3) round(3/5*scrsz(4))]); % Full screen figure
    for id_axe=1:1:3
        sub_axes=subplot(1,3,id_axe,'Parent',Fig_);
        hold(sub_axes,'on');
        title (sprintf('Averaged grey level, view normal to axe %i',id_axe),'FontName','Times New Roman','FontSize',16);
        if id_axe==1
            [Y,X]=meshgrid(1:1:domain_size(2),1:1:domain_size(3)); % Create grid
            xlabel('3rd Axis'); ylabel('2nd Axis');
        elseif id_axe==2
            [Y,X]=meshgrid(1:1:domain_size(1),1:1:domain_size(3));
            xlabel('3rd Axis'); ylabel('1st Axis');
        elseif id_axe==3
            [Y,X]=meshgrid(1:1:domain_size(1),1:1:domain_size(2));
            xlabel('2nd Axis'); ylabel('1st Axis');
        end
        surfc(X,Y,greylevel_normal(id_axe).map','Parent',sub_axes,'LineStyle','none'); % Create surface plot with contour
        colorbar('peer',sub_axes); % Create colorbar
        set(sub_axes,'FontName','Times New Roman','FontSize',14);
        axis(sub_axes,'equal');
        axis(sub_axes,'tight');
        % View
        view(0,90)
    end

else
    Fig_=[];
end

end

