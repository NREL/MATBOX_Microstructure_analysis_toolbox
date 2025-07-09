function [outputImage,newtype,numLabels] = fct_superpixels(M,p)

sz = size(M);
dimension = length(sz);

if strcmp(p.NumChoice,'superpixels')
    p.numlabel = p.num;
elseif strcmp(p.NumChoice,'pixels per superpixel (volume)')
    p.numlabel = round( numel(M)/p.num );
elseif strcmp(p.NumChoice,'pixels per superpixel (length)')
    pixel_per_superpixel = p.num ^ dimension;
    p.numlabel = round( numel(M)/pixel_per_superpixel );
end

if p.Refine
    p.Method = "slic0";
else
    p.Method = "slic";
end

newtype = 'same';

if dimension == 2
    if strcmp(p.Compactness_InitialValue_choice,'Default')
        [L,numLabels] = superpixels(M, p.numlabel, NumIterations=p.NumIterations, Method=p.Method);
    elseif strcmp(p.Compactness_InitialValue_choice,'Custom')
        [L,numLabels] = superpixels(M, p.numlabel, NumIterations=p.NumIterations, Method=p.Method, Compactness=p.Compactness_InitialValue);
    end

    % Assign each superpixel to the mean
    outputImage = zeros(size(M),'like',M);
    idx = label2idx(L);
    for labelVal = 1:numLabels
        Idx = idx{labelVal};
        outputImage(Idx) = mean(M(Idx));
    end

else
    if strcmp(p.Compactness_InitialValue_choice,'Default')
        [L,numLabels] = superpixels3(M, p.numlabel, NumIterations=p.NumIterations, Method=p.Method);
    elseif strcmp(p.Compactness_InitialValue_choice,'Custom')
        [L,numLabels] = superpixels3(M, p.numlabel, NumIterations=p.NumIterations, Method=p.Method, Compactness=p.Compactness_InitialValue);
    end

    % Assign each superpixel to the mean
    pixelIdxList = label2idx(L);
    outputImage = zeros(size(M),'like',M);
    for superpixel = 1:numLabels
        memberPixelIdx = pixelIdxList{superpixel};
        outputImage(memberPixelIdx) = mean(M(memberPixelIdx));
    end

end

[outputImage] = fct_intconvert(outputImage);

end