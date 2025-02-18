function [M,newtype,foo] = fct_texture(M,p)

sz = size(M);
foo = [];
newtype = 'Channel';
dimension = length(sz);

if dimension==2
    if strcmp(p.Morphologicalstructuringelement,'Square or cube')
        SE = strel("square",2*p.range+1);
    elseif strcmp(p.Morphologicalstructuringelement,'Disc or sphere')
        SE = strel("disk",2*p.range+1);
    end
else
    if strcmp(p.Morphologicalstructuringelement,'Square or cube')
        SE = strel("cube",2*p.range+1);
    elseif strcmp(p.Morphologicalstructuringelement,'Disc or sphere')
        SE = strel("sphere",2*p.range+1);
    end
end

if strcmp(p.metric,'Standard deviation')
    M = stdfilt(M,SE.Neighborhood);
elseif strcmp(p.metric,'Entropy')
    M = entropyfilt(M,SE.Neighborhood);
elseif strcmp(p.metric,'Max-min range')
    M = rangefilt(M,SE.Neighborhood);
end

%[M] = fct_intconvert(M);

end
