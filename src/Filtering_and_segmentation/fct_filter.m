function [M,newtype,foo] = fct_filter(M,p)

sz = size(M);
foo=[];
newtype = 'same';
dimension = length(sz);

if strcmp(p.filter,'Mean filter')
    if dimension==2
        M = imboxfilt(M,p.filtersize,'Padding',p.padding);
    else
        M = imboxfilt3(M,p.filtersize,'Padding',p.padding);
    end
elseif strcmp(p.filter,'Gaussian filter')
    if dimension==2
        M = imgaussfilt(M,p.sigma,'FilterSize',p.filtersize,'Padding',p.padding,'FilterDomain',p.domain);
    else
        M = imgaussfilt3(M,p.sigma,'FilterSize',p.filtersize,'Padding',p.padding,'FilterDomain',p.domain);
    end
end

M = round(M);

% Convert to uint8 or uint16
[M] = fct_intconvert(M);

end