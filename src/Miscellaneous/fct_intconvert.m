function [M] = fct_intconvert(M)

decimalCheck = mod(M, 1) ~= 0;
if sum(sum(sum(decimalCheck))) >=1
    M = double(M);
    min_ = min(min(min( M )));
    max_ = max(max(max( M )));
    new_max = 2^16-1;
    a = (new_max-0)/(max_-min_);
    b = 0-a*min_;
    M = round(a.*M + b);
end

if min(min(min(M)))>=0
    if max(max(max(M)))<=2^8-1
        M=uint8(M);
    elseif max(max(max(M)))<=2^16-1
        M=uint16(M);
    elseif max(max(max(M)))<=2^32-1
        M=uint32(M);        
    end
else
    if min(min(min(abs(M))))>= - 2^7 && max(max(max(abs(M)))) <= 2^7 -1 
        M=int8(M);
    elseif min(min(min(abs(M))))>= - 2^15 && max(max(max(abs(M)))) <= 2^15 -1 
        M=int16(M);
    elseif min(min(min(abs(M))))>= - 2^31 && max(max(max(abs(M)))) <= 2^31 -1 
        M=int32(M);        
    end
end

end