function [M,newtype,p] = fct_Crop(M,p)
sz = size(M);
dimension = length(sz);
newtype = 'same';

if p.autoset_tomatch_downscaling
    if p.downscaling_factor>1
        cropxyz = zeros(dimension,1);
        for k=1:1:dimension
            tmp = [1:1:sz(k)];
            tmp2 = mod(tmp,p.downscaling_factor);
            idx = max(find(tmp2==0));
            cropxyz(k) = length(tmp)-idx;
        end
        p.x0 = ones(dimension,1);
        p.x1 = sz'-cropxyz;
        p.rel = false;
    end
end

if p.rel
    x0sav = p.x0;
    x1sav = p.x1;
    a = (sz-1)./(1-0);
    b = 1;
    p.x0 = round(a'.*p.x0 + b');
    p.x1 = round(a'.*p.x1 + b');
end

if p.DeltaROI
    p.x1 = sz' - p.delta_end;
end

if p.x1(1)>sz(1) || p.x1(2)>sz(2)
    f = warndlg('Cannot crop! Volume is too small','Warning');
    return
end

if dimension == 2
    M = M(p.x0(1):p.x1(1), p.x0(2):p.x1(2));
else
    if p.x1(3)>sz(3)
        f = warndlg('Cannot crop! Volume is too small','Warning');
        return
    end
    M = M(p.x0(1):p.x1(1), p.x0(2):p.x1(2), p.x0(3):p.x1(3));
end

sz = size(M);
if length(sz)>2 && sum(sz==1) && p.autoreshape
    sz(sz==1)=[];
    M=reshape(M,sz);
end

if p.rel
    p.x0 = x0sav;
    p.x1 = x1sav;
end

end