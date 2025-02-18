function [M,newtype,foo] = fct_highlow_filter(M,p)

sz = size(M);
foo=[];
newtype = 'same';

if strcmp(p.choice,'equal or above')
    M(M>=p.threshold) = p.newval;
elseif strcmp(p.choice,'equal or below')
    M(M<=p.threshold) = p.newval;
end

[M] = fct_intconvert(M);

end