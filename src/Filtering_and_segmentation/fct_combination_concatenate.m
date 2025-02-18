function [M,newtype] = fct_combination_concatenate(Ma,Mb,p)

success = false;

if strcmp(p.Ma_type,p.Mb_type)
    newtype = p.Ma_type;
elseif (strcmp(p.Ma_type,'Segmented (phase)') && strcmp(p.Mb_type,'Segmented (instance)')) || (strcmp(p.Ma_type,'Segmented (instance)') && strcmp(p.Mb_type,'Segmented (phase)'))
    newtype = 'Segmented (instance)';
else
    newtype = 'Grey level';   
end

szA = size(Ma);
dimA = length(szA);

szB = size(Mb);
dimB = length(szB);

Ma = double(Ma);
Mb = double(Mb);

if p.axe == 1
    if dimA == dimB
        M=[Ma;Mb];
        success = true;
    end

elseif p.axe == 2
    if dimA == dimB
        M=[Ma Mb];
        success = true;
    end

else
    if dimA==2
        szA(3) = 1;
    end
    if dimB==2
        szB(3) = 1;
    end
    if szA(1)==szB(1) && szA(2)==szB(2)
        M = zeros(szA(1),szA(2),szA(3)+szB(3));
        M(:,:,1:szA(3)) = Ma;
        M(:,:,szA(3)+1:end) = Mb;
        success = true;
    end
end
  
if success
    [M] = fct_intconvert(M);
else
    M = 0;
end

end