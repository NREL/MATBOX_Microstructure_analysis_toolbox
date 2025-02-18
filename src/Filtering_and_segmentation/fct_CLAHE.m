function [M,newtype,foo] = fct_CLAHE(M,p)

sz = size(M);
foo=[];
newtype = 'same';
dimension = length(sz);
if dimension==2
    sz = [sz 1];
end

for z=1:1:sz(3)
    if strcmp(p.Distribution,'uniform')
        M(:,:,z) = adapthisteq(M(:,:,z),'NumTiles',p.NumTiles,'ClipLimit',p.ClipLimit,'NBins',p.NBins,'Range',p.Range,'Distribution',p.Distribution);
    else
        M(:,:,z) = adapthisteq(M(:,:,z),'NumTiles',p.NumTiles,'ClipLimit',p.ClipLimit,'NBins',p.NBins,'Range',p.Range,'Distribution',p.Distribution,'Alpha',p.Alpha);
    end
end

% Convert to uint8 or uint16
[M] = fct_intconvert(M);

end