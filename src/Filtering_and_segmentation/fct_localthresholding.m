function [Mseg,newtype,allthresholds] = fct_localthresholding(Mgrey,p)

newtype = 'Segmented (phase)';

if strcmp(p.option,'Otsu')
    [Mseg,allthresholds] = fct_localthresholding_Otsu(Mgrey,p);
elseif strcmp(p.option,'Manual')
    [Mseg,allthresholds] = fct_localthresholding_Manual(Mgrey,p);
end

end