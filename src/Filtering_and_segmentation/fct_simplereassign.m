function [M,newtype,foo] = fct_simplereassign(M,p)
% This function only exist as it allows use with feval in macros

foo=[];
if strcmp(p.action,'Absolute')
    M = abs(M);
    newtype = 'same';
elseif strcmp(p.action,'Round')
    M = round(M,p.val);
    newtype = 'same';
elseif strcmp(p.action,'Arithmetic')
    if strcmp(p.ope,'+')
        M = double(M) + double(p.val);    
    elseif strcmp(p.ope,'x')
        M = double(M) .* double(p.val);
    end
    newtype = 'same';
    [M] = fct_intconvert(M);
elseif strcmp(p.action,'Not')
    M=~M;      
    newtype = 'Segmented (phase)';
    [M] = fct_intconvert(M);
elseif strcmp(p.action,'Binarize')
    %M(M~=0)=1;  
    tmp = M;
    newtype = 'Segmented (phase)';  
    if strcmp(p.choice,'M( M ~= x) = 1, 0 otherwise')
        
        M = zeros(size(M));
        M(tmp~=p.val) = 1;
    elseif strcmp(p.choice,'M( M == x) = 1, 0 otherwise')
        M = zeros(size(M));
        M(tmp==p.val) = 1;
    elseif strcmp(p.choice,'M( M ~= x) = 0, 1 otherwise')
        M = ones(size(M));
        M(tmp~=p.val) = 0;
    elseif strcmp(p.choice,'M( M == x) = 0, 1 otherwise')
        M = ones(size(M));
        M(tmp==p.val) = 0;        
    elseif strcmp(p.choice,'M( M ==min(M) ) = 1, 0 otherwise') 
        M = zeros(size(M));
        M( tmp==min(min(min(tmp))) ) = 1;
    elseif strcmp(p.choice,'M( M ==max(M) ) = 1, 0 otherwise')  
        M = zeros(size(M));
        M( tmp==max(max(max(tmp))) ) = 1;   
    end
    
    [M] = fct_intconvert(M);


elseif strcmp(p.action,'Substitution')
    newtype = 'same';  
    if strcmp(p.choice,'M( M == x) = y')
        M( M==p.x ) = p.y;
    elseif strcmp(p.choice,'M( M ~= x) = y')
        M( M~=p.x ) = p.y;        
    elseif strcmp(p.choice,'M( M == min(M) ) = y')
        M( M==min(min(min(M))) ) = p.y;
    elseif strcmp(p.choice,'M( M == max(M) ) = y')
        M( M==max(max(max(M))) ) = p.y;   
    elseif strcmp(p.choice,'M( M > x) = y')
        M( M>p.x ) = p.y;   
    elseif strcmp(p.choice,'M( M < x) = y')
        M( M<p.x ) = p.y;  
    end    
    [M] = fct_intconvert(M);
end

end