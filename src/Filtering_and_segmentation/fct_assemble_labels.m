function [M,newtype,foo] = fct_assemble_labels(M,p)

newtype = 'Segmented (instance)';
foo = 1;

sz = size(M);
dimension = length(sz);

n_parts = length(p.parts);
x_min = 1e9; x_max = -1e9;
y_min = 1e9; y_max = -1e9;
if dimension == 3
    z_min = 1e9; z_max = -1e9;
else
    z_min = 1; z_max = 1;
end
for k = 1:n_parts
    idx = find(M==p.parts(k));
    if dimension == 2
        [IX,IY] = ind2sub(sz,idx);
    else
        [IX, IY, IZ] = ind2sub(sz,idx);
    end
    x_min = min([x_min min(IX)]);
    y_min = min([y_min min(IY)]);
    x_max = max([x_max max(IX)]);
    y_max = max([y_max max(IY)]);
    if dimension == 3
        z_min = min([z_min min(IZ)]);
        z_max = max([z_max max(IZ)]);
    end
end

sub = M(x_min:x_max, y_min:y_max, z_min:z_max);
sz_sub = size(sub);
BWn = zeros(sz_sub);
for k = 1:n_parts
    BWn(sub==p.parts(k))=k;
    sub(sub==p.parts(k))=0;
end

% Fill voids
BW = BWn;
BW(BW~=0)=1;
if strcmp(p.method,'Convex hull')
    if dimension==3
        res = regionprops3(BW,'ConvexImage');
        res = cell2mat(res.ConvexImage);
    else
        res = regionprops(BW(:,:,5),'ConvexImage');
        res = double(res.ConvexImage);
    end
elseif strcmp(p.method,'Sphere-based local dilation')
    pp.choice = 'One label';
    pp.labeltodilate = 1;
    pp.removeotherlabel = true;
    pp.background = 0;
    pp.per_cluster = false;
    pp.per_label = false;
    pp.spherediameter = p.diameterthreshold;
    pp.round_dmap = true;
    pp.extendFOV = true;
    res = fct_localdilation_spheresizefilter_algo(BW,pp);
end

if p.overwritte
    sub(res==1) = p.parts(1);
else
    cond1 = res==1;
    cond2 = sub==0;
    sub(cond1+cond2==2) = p.parts(1);
end

% Replace
M(x_min:x_max, y_min:y_max, z_min:z_max) = sub;

end