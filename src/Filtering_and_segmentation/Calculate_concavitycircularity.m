function [Asolidity, Acircularity] = Calculate_concavitycircularity(BW,Edges,p)

sz = size(Edges);
dimension = length(sz);

if strcmp(p.metric,'solidity')
    Asolidity = zeros(sz);
else
    Asolidity = [];
end
if strcmp(p.metric,'circularity')
    Acircularity = zeros(sz);
else
    Acircularity = [];
end

idx = find(Edges==2);
n = length(idx);
if dimension==2
    [IX, IY] = ind2sub(sz,idx);
    x_mins = max([ones(n,1) IX-p.wndrange],[],2);
    y_mins = max([ones(n,1) IY-p.wndrange],[],2);
    x_maxs = min([ones(n,1)*sz(1) IX+p.wndrange],[],2);
    y_maxs = min([ones(n,1)*sz(2) IY+p.wndrange],[],2);

    for k=1:1:n
        BWsub = BW(x_mins(k):x_maxs(k), y_mins(k):y_maxs(k));
        C = bwlabel(BWsub);
        if length(unique(C))>2
            if IX(k)-p.wndrange>=1 && IX(k)+p.wndrange<=sz(1)
                x = p.wndrange+1;
            elseif IX(k)-p.wndrange<1 && IX(k)+p.wndrange<=sz(1)
                x = IX(k);
            elseif IX(k)-p.wndrange>=1 && IX(k)+p.wndrange>sz(1)
                x = p.wndrange+1;
            elseif IX(k)-p.wndrange<1 && IX(k)+p.wndrange>sz(1)
                x = IX(k);
            end

            if IY(k)-p.wndrange>=1 && IY(k)+p.wndrange<=sz(2)
                y = p.wndrange+1;
            elseif IY(k)-p.wndrange<1 && IY(k)+p.wndrange<=sz(2)
                y = IY(k);
            elseif IY(k)-p.wndrange>=1 && IY(k)+p.wndrange>sz(2)
                y = p.wndrange+1;
            elseif IY(k)-p.wndrange<1 && IY(k)+p.wndrange>sz(2)
                y = IY(k);
            end 
            BWsub(C~=C(x,y))=0;  
        end

        if strcmp(p.metric,'solidity')
            Solidity = regionprops(BWsub,"Solidity");
            Solidity = min([Solidity.Solidity, 1]);
            Asolidity(IX(k), IY(k)) = 1-Solidity;
        end
        if strcmp(p.metric,'circularity')
            Circularity = regionprops(BWsub,"Circularity");
            Circularity = min([Circularity.Circularity, 1]);
            Acircularity(IX(k), IY(k)) = 1-Circularity;
        end

        %statsConvexHull = regionprops(BWsub,"ConvexHull")
        %statsConvexImage = regionprops(BWsub,"ConvexImage");

        % statsConvexArea = regionprops(BWsub,"ConvexArea");
        % nsub = sum(sum(BWsub));
        % if p.norm_local
        %     Concavity(IX(k), IY(k)) = (statsConvexArea.ConvexArea - nsub)/nsub;
        % else
        %     Concavity(IX(k), IY(k)) = statsConvexArea.ConvexArea - nsub;
        % end
    end

else
    [IX, IY, IZ] = ind2sub(sz,idx);
    x_mins = max([ones(n,1) IX-p.wndrange],[],2);
    y_mins = max([ones(n,1) IY-p.wndrange],[],2);
    z_mins = max([ones(n,1) IZ-p.wndrange],[],2);
    
    x_maxs = min([ones(n,1)*sz(1) IX+p.wndrange],[],2);
    y_maxs = min([ones(n,1)*sz(2) IY+p.wndrange],[],2);
    z_maxs = min([ones(n,1)*sz(3) IZ+p.wndrange],[],2);
    
    for k=1:1:n
        BWsub = BW(x_mins(k):x_maxs(k), y_mins(k):y_maxs(k), z_mins(k):z_maxs(k));
        C = bwlabeln(BWsub);
        if length(unique(C))>2
            if IX(k)-p.wndrange>=1 && IX(k)+p.wndrange<=sz(1)
                x = p.wndrange+1;
            elseif IX(k)-p.wndrange<1 && IX(k)+p.wndrange<=sz(1)
                x = IX(k);
            elseif IX(k)-p.wndrange>=1 && IX(k)+p.wndrange>sz(1)
                x = p.wndrange+1;
            elseif IX(k)-p.wndrange<1 && IX(k)+p.wndrange>sz(1)
                x = IX(k);
            end

            if IY(k)-p.wndrange>=1 && IY(k)+p.wndrange<=sz(2)
                y = p.wndrange+1;
            elseif IY(k)-p.wndrange<1 && IY(k)+p.wndrange<=sz(2)
                y = IY(k);
            elseif IY(k)-p.wndrange>=1 && IY(k)+p.wndrange>sz(2)
                y = p.wndrange+1;
            elseif IY(k)-p.wndrange<1 && IY(k)+p.wndrange>sz(2)
                y = IY(k);
            end 

            if IZ(k)-p.wndrange>=1 && IZ(k)+p.wndrange<=sz(3)
                z = p.wndrange+1;
            elseif IZ(k)-p.wndrange<1 && IZ(k)+p.wndrange<=sz(3)
                z = IZ(k);
            elseif IZ(k)-p.wndrange>=1 && IZ(k)+p.wndrange>sz(3)
                z = p.wndrange+1;
            elseif IZ(k)-p.wndrange<1 && IZ(k)+p.wndrange>sz(3)
                z = IZ(k);
            end 

            BWsub(C~=C(x,y,z))=0;  
        end

        if strcmp(p.metric,'solidity')
            Solidity = min([regionprops3(BWsub,"Solidity"), 1]);
            Asolidity(IX(k), IY(k), IZ(k)) = 1-Solidity.Solidity;
        end
        if strcmp(p.metric,'circularity')
            % statsCircularity = regionprops(BWsub,"Circularity"); 2D only
            V = sum(sum(sum(BWsub)));
            S = regionprops3(BWsub,"SurfaceArea");
            S = S.SurfaceArea;
            %[Edges] = Find_edgesBW(BWsub);
            %S = sum(sum(sum(Edges==2)));
            % sphericity = 9/(2*S) * V^(2/3) * (4*pi/3)^(1/3);
            %sphericity = pi^(1/3) * (6*V)^(2/3) / S; % Values are too high
            %sphericity = pi^(1/3) * (6*V)^(2/3) / (2/3*S) * 0.5;
            sphericity = min([pi^(1/3) * (6*V)^(2/3) / S, 1]);
            Acircularity(IX(k), IY(k), IZ(k)) = 1-sphericity;
        end

        %statsConvexHull = regionprops(BWsub,"ConvexHull")
        %statsConvexImage = regionprops(BWsub,"ConvexImage");

        % statsConvexArea = regionprops(BWsub,"ConvexArea");
        % nsub = sum(sum(BWsub));
        % if p.norm_local
        %     Concavity(IX(k), IY(k)) = (statsConvexArea.ConvexArea - nsub)/nsub;
        % else
        %     Concavity(IX(k), IY(k)) = statsConvexArea.ConvexArea - nsub;
        % end
    end


end

% if ~p.norm_local
%     Concavity = Concavity/max(max(max(Concavity)));
% end

% Round and filter to smooth and remove numerical error induced by pixel representation
if strcmp(p.metric,'solidity')
    Asolidity = round(Asolidity,2);
    Asolidity(Asolidity<=p.threshold_metric)=0;
end
if strcmp(p.metric,'circularity')
    Acircularity = round(Acircularity,2);
    Acircularity(Acircularity<=p.threshold_metric)=0;
end

end