function [M,newtype,foo] = fct_linearrescaling(M,p)

sz = size(M);
foo=[];
newtype = 'same';

max_=max(max(max(M)));
min_=min(min(min(M)));
initial_delta=max_-min_;
final_delta=p.nvalues-1;
M =  round( ((double((M-min_)) ./ double(initial_delta)) .* final_delta)+1 );

% Convert to uint8 or uint16
[M] = fct_intconvert(M);

end