function [M,newtype] = fct_combination_difference(Ma,Mb,p)

if strcmp(p.Ma_type,p.Mb_type)
    newtype = p.Ma_type;
elseif (strcmp(p.Ma_type,'Segmented (phase)') && strcmp(p.Mb_type,'Segmented (instance)')) || (strcmp(p.Ma_type,'Segmented (instance)') && strcmp(p.Mb_type,'Segmented (phase)'))
    newtype = 'Segmented (instance)';
else
    newtype = 'Grey level';   
end

if strcmp(p.difference,'A-B')
    M = Ma - Mb;
elseif strcmp(p.difference,'abs(A-B)')
    M = abs(Ma - Mb);
elseif strcmp(p.difference,'(A-B) / ((A+B)/2)')
    M = (Ma-Mb) ./ ((Ma+Mb)/2);
elseif strcmp(p.difference,'abs( (A-B) / ((A+B)/2) )')
    M = abs((Ma-Mb) ./ ((Ma+Mb)/2));
elseif strcmp(p.difference,'(A-B) / <(A+B)/2)>')
    deno = sum(sum(sum( Ma+Mb ))) / (2*numel(Ma));
    M = (Ma-Mb) / deno;
elseif strcmp(p.difference,'abs( (A-B) / <(A+B)/2)> )')
    deno = sum(sum(sum( Ma+Mb ))) / (2*numel(Ma));
    M = abs((Ma-Mb) / deno);
end

if sum(sum(sum( M==round(M) )))==numel(M)
    [M] = fct_intconvert(M);
end

end