function [M,newtype] = fct_combination_multiplication(Mseg,Mlabel,p)

checkmissingvoxels_is_relevant = false;
if strcmp(p.Ma_type,p.Mb_type)
    newtype = p.Ma_type;
elseif (strcmp(p.Ma_type,'Segmented (phase)') && strcmp(p.Mb_type,'Segmented (instance)')) || (strcmp(p.Ma_type,'Segmented (instance)') && strcmp(p.Mb_type,'Segmented (phase)'))
    newtype = 'Segmented (instance)';
    checkmissingvoxels_is_relevant = true;
else
    newtype = 'Grey level';   
end

M = Mseg .* Mlabel;

if p.checkmissingvoxels && checkmissingvoxels_is_relevant % Missing points ? Assigned to nearest labels
    if strcmp(p.Ma_type,'Segmented (instance)')
        tmp = Mseg;
        Mseg = Mlabel;
        Mlabel = tmp;
    end
    sz = size(M);
    tmp = zeros(sz);
    tmp(Mlabel==0)=1;
    missingpoints = Mseg.*double(tmp);
    id_missingpoints = find(missingpoints);
    if ~isempty(id_missingpoints)>0 % Yes, assign to nearest label lake
        [~,idx] = bwdist(~tmp);
        M(id_missingpoints) = Mlabel(idx(id_missingpoints));
    end
end

[M] = fct_intconvert(M);

end

