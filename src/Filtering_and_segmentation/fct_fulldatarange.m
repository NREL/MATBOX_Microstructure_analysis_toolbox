function [M,newtype,foo] = fct_fulldatarange(M,p)

sz = size(M);
foo=[];
newtype = 'same';

min_=min(min(min(M)));
M = M - min_;

max_=max(max(max(M)));

if max_<=2^8-1
    up = 2^8-1;
else
    up = 2^16-1;
end
a = up/max_;

M = round(a*M);

% Format
[M] = fct_intconvert(M);

end