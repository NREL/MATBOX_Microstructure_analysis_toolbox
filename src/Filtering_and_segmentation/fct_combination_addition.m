function [M,newtype] = fct_combination_addition(Ma,Mb,p)

if strcmp(p.Ma_type,p.Mb_type)
    newtype = p.Ma_type;
elseif (strcmp(p.Ma_type,'Segmented (phase)') && strcmp(p.Mb_type,'Segmented (instance)')) || (strcmp(p.Ma_type,'Segmented (instance)') && strcmp(p.Mb_type,'Segmented (phase)'))
    newtype = 'Segmented (instance)';
else
    newtype = 'Grey level';   
end

if p.plus
    M = double(Ma) + max(max(max(double(Mb))));
else
    M = double(Ma) + double(Mb);
end

if sum(sum(sum( M==round(M) )))==numel(M)
    [M] = fct_intconvert(M);
end

end