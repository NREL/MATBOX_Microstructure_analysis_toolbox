function [M2,newtype,foo] = fct_valfilter(M,p)

sz = size(M);
foo=[];
newtype = 'Segmented (phase)';
dimension = length(sz);

if p.scalewithdim
    p.threshold = p.threshold ^ dimension;
end

M2 = zeros(sz);
if strcmp(p.choice,'<=')
    M2(M <= p.threshold) = 1;    
elseif strcmp(p.choice,'>')
    M2(M > p.threshold) = 1;   
end

if p.background>=0
    M2(M==p.background)=0;
end

M2 = round(M2);

% Convert to uint8 or uint16
[M2] = fct_intconvert(M2);

end