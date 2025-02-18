function [M,newtype,foo] = fct_histequalization(M,p)

sz = size(M);
foo=[];
newtype = 'same';
dimension = length(sz);

% Convert to uint8 or uint16 (histeq does not work with double)
if min(min(min(M)))>=0
    if max(max(max(M)))<=255
        M=uint8(M);
    elseif max(max(max(M)))<=65535
        M=uint16(M);
    end
else
    if max(max(max(abs(M))))<=255
        M=int8(M);
    elseif max(max(max(abs(M))))<=65535
        M=int16(M);
    end
end

if ~isempty(p.nbin)
    if p.relative_nbin
        p.nbin = round( p.ratio_nbin * length(unique(M)) );
    end
    if p.perslice && dimension==3
        for z=1:1:sz(3)
            M(:,:,z) = histeq(M(:,:,z),p.nbin);
        end
    else
        M = histeq(M,p.nbin);
    end

elseif ~isempty(p.hgram)
    if p.perslice && dimension==3
        for z=1:1:sz(3)
            M(:,:,z) = histeq(M(:,:,z),p.hgram);
        end
    else
        M = histeq(M,p.hgram);
    end

elseif ~isempty(p.M2)
    if p.scalerange
        min_ = min(min(min(M)));
        max_ = max(max(max(M)));
        n = max_-min_+1;
        c = zeros(n,1);
        for k=1:1:n
            c(k) = sum(sum(sum(M==min_+k-1)));
        end
        M = histeq(M,c);
    else % Data range
        if min(min(min(M)))>=0 && max(max(max(M)))<=255
            min_ = 0;
            max_ = 255;
        elseif min(min(min(M)))>=0 && max(max(max(M)))<=65535
            min_ = 0;
            max_ = 65535;
        end
        n = max_-min_+1;
        c = zeros(n,1);
        for k=1:1:n
            c(k) = sum(sum(sum(M==min_+k-1)));
        end
        M = histeq(M,c);
    end
end


% Convert to uint8 or uint16
if min(min(min(M)))>=0
    if max(max(max(M)))<=255
        M=uint8(M);
    elseif max(max(max(M)))<=65535
        M=uint16(M);
    end
else
    if max(max(max(abs(M))))<=255
        M=int8(M);
    elseif max(max(max(abs(M))))<=65535
        M=int16(M);
    end
end

end