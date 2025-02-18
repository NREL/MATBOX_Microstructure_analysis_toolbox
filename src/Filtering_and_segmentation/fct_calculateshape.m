function [Ms,newtype,foo] = fct_calculateshape(M,p)

newtype = 'Channel';
foo = 1;

sz = size(M);
dimension = length(sz);
Ms = zeros(sz,'single');

M = M+1; % no zero for regionprops

if strcmp(p.choice,'All labels')
    labels = unique(M);
    n = max(max(max(M)));
elseif strcmp(p.choice,'One label')
    tmp = zeros(sz);
    tmp(M==p.label + 1)=1;
    M=tmp;
    clear tmp;
    labels = 1;
    n = 1;
elseif strcmp(p.choice,'All labels except some labels')
    labels = unique(M);
    n = max(max(max(M)));
    excluded_labels = str2num(p.excluded_labels);
    if ~isempty(excluded_labels)
        excluded_labels = excluded_labels + 1;
        for k=1:1:length(excluded_labels)
            labels(labels==excluded_labels(k))=[];
        end
    end
end

unit = ones(n,1);
res = zeros(n,1);


if dimension == 2
    idx_tmp = regionprops(M,"PixelIdxList");

    if strcmp(p.metric,'Circularity (2D), sphericity (3D)')
        res = regionprops("table",M,"Circularity");
        res = min([res.Circularity, unit],[],2);
 
    % elseif strcmp(p.metric,'Eccentricity')    
    %     res = regionprops("table",M,"Eccentricity");
    %     res = res.Eccentricity;
    % 
    % elseif strcmp(p.metric,'EulerNumber')    
    %     res = regionprops("table",M,"EulerNumber");
    %     res = res.EulerNumber;   
    % 
    %  elseif strcmp(p.metric,'Orientation')    
    %     res = regionprops("table",M,"Orientation");
    %     res = res.Orientation;          

    elseif strcmp(p.metric,'Solidity')
        res = regionprops("table",M,"Solidity");
        res = min([res.Solidity, unit],[],2);   

    elseif strcmp(p.metric,'Aspect ratio')
        k2 = 0;
        for k=1:1:n
            if sum(labels==k)
                k2=k2+1;
                idx = idx_tmp.VoxelIdxList{labels(k2),1};
                [IX,IY] = ind2sub(sz,idx);
                Lx = max(IX) - min(IX) + 1;
                Ly = max(IY) - min(IY) + 1;
                res(k) = min([Lx Ly])/max([Lx Ly]);
            end
        end       

    end
else
    idx_tmp = regionprops3(M,"VoxelIdxList");

    if strcmp(p.metric,'Circularity (2D), sphericity (3D)')
        V = regionprops3(M,"Volume");
        V = V.Volume;
        S = regionprops3(M,"SurfaceArea");
        S = S.SurfaceArea;
        res = min([(pi^(1/3) .* (6*V).^(2/3)) ./ S, unit],[],2);
        
    elseif strcmp(p.metric,'Solidity')
        % res = regionprops3("table",M,"Solidity");
        res = regionprops3(M,"Solidity");
        res = min([res.Solidity, unit],[],2);

    elseif strcmp(p.metric,'Aspect ratio')
        k2 = 0;
        for k=1:1:n
            if sum(labels==k)
                k2=k2+1;
                idx = idx_tmp.VoxelIdxList{labels(k2),1};
                [IX,IY,IZ] = ind2sub(sz,idx);
                Lx = max(IX) - min(IX) + 1;
                Ly = max(IY) - min(IY) + 1;
                Lz = max(IZ) - min(IZ) + 1;
                res(k) = min([Lx Ly Lz])/max([Lx Ly Lz]);
            end
        end

    end
end


n2 = length(labels);
if dimension == 2
    for k = 1:1:n2
        idx = idx_tmp(labels(k)).PixelIdxList;
        Ms(idx) = res(labels(k));
    end
else
    for k = 1:1:n2
        idx = idx_tmp.VoxelIdxList{labels(k),1};
        Ms(idx) = res(labels(k));
    end
end

if p.convert_toint
    Ms = round(Ms*power(10,p.precision));
    % Convert to uint8 or uint16
    [Ms] = fct_intconvert(Ms);
end

end