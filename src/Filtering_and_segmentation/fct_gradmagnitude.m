function [Gmag,newtype,p] = fct_gradmagnitude(M,p)

sz = size(M);
newtype = 'Channel';
dimension = length(sz);

if dimension==2
    [Gmag,~] = imgradient(M,p.method);
else
    if strcmp(p.method,'roberts') % Only available for 2D
        p.method = 'sobel';
    end
    [Gmag,~,~] = imgradient3(M,p.method);
end

end
