function [Ms,newtype,foo] = fct_calculateconnectivity(M,p)

newtype = 'Segmented (phase)';
foo = 1;

sz = size(M);
dimension = length(sz);

Ms = zeros(sz);

if dimension == 2
    conn_shape = zeros(2,2);
    ncheck = 4;
    p.conditions = p.conditions(:,1:2);
else
    conn_shape = zeros(2,3);
    ncheck = 6;
end

if strcmp(p.mustverify,'All conditions')
    if p.onlycheckpositive
        ntarget = sum(sum(p.conditions));
    else
        ntarget = ncheck;
    end
elseif strcmp(p.mustverify,'At least one condition')
    ntarget = 1;
end

labels = unique(M);
labels(labels==p.background)=[];
n = length(labels);

for k = 1:1:n
    conn = conn_shape;
    idx = find(M==labels(k));
    if dimension == 2
        [IX,IY] = ind2sub(sz,idx);
    else
        [IX,IY,IZ] = ind2sub(sz,idx);
        if min(IZ)==1
            conn(1,3)=1;
        end
        if max(IZ)==sz(3)
            conn(2,3)=1;
        end
    end
    if min(IX)==1
        conn(1,1)=1;
    end
    if min(IY)==1
        conn(1,2)=1;
    end
    if max(IX)==sz(1)
        conn(2,1)=1;
    end
    if max(IY)==sz(2)
        conn(2,2)=1;
    end

    idequal = find(conn==p.conditions);
    ntot = length(idequal);
    nplus = sum(sum(conn(idequal)));
    if p.onlycheckpositive
        if nplus>=ntarget
            Ms(idx) = 1;
        end
    else
        if ntot>=ntarget
            Ms(idx) = 1;
        end
    end    
end

if p.inverse
    Ms = double(~Ms);
end

[Ms] = fct_intconvert(Ms);

end