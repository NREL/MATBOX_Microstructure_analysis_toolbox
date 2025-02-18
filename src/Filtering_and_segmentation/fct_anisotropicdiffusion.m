function [M,newtype,str_p] = fct_anisotropicdiffusion(M,p)
sz = size(M);
dimension = length(sz);
newtype = 'same';

if dimension==2
    sz = [sz 1];
end

if strcmp(p.apply,'Single slice') || dimension==2
    if p.estimategradientiteration
        [Gradient_threshold,number_iteration] = imdiffuseest(M,'Connectivity',p.connectivity,'ConductionMethod',p.conduction);
        str_p = ['Number of iterations: ' num2str(number_iteration,'%i') ', gradient thershold for each iteration: ' num2str(Gradient_threshold)];
        M = imdiffusefilt(M,'Connectivity',p.connectivity,'ConductionMethod',p.conduction,'GradientThreshold',Gradient_threshold,'NumberOfIterations',number_iteration);
    else
        Gradient_threshold = ones(1,p.number_iterations)*p.gradthreshold;
        str_p = ['Number of iterations: ' num2str(p.number_iterations,'%i') ', gradient thershold for each iteration: ' num2str(Gradient_threshold)];
        M = imdiffusefilt(M,'Connectivity',p.connectivity,'ConductionMethod',p.conduction,'GradientThreshold',Gradient_threshold,'NumberOfIterations',p.number_iterations);
    end

elseif strcmp(p.apply,'Full volume')
    Gradient_threshold = ones(1,p.number_iterations)*p.gradthreshold;
    str_p = ['Applied full volume at once. Number of iterations: ' num2str(p.number_iterations,'%i') ', gradient thershold for each iteration: ' num2str(Gradient_threshold)];
    M = imdiffusefilt(M,'Connectivity',p.connectivity,'ConductionMethod',p.conduction,'GradientThreshold',Gradient_threshold,'NumberOfIterations',p.number_iterations);

elseif strcmp(p.apply,'Slice per slice, along axis 1')
    for x=1:1:sz(1)
        slice_ = squeeze(M(x,:,:));
        if p.estimategradientiteration
            str_p = 'Applied per slice along axis 1, auto estimation of number of iterations and gradient threshold';
            [Gradient_threshold,number_iteration] = imdiffuseest(slice_,'Connectivity',p.connectivity,'ConductionMethod',p.conduction);
            slice_ = imdiffusefilt(slice_,'Connectivity',p.connectivity,'ConductionMethod',p.conduction,'GradientThreshold',Gradient_threshold,'NumberOfIterations',number_iteration);
        else
            Gradient_threshold = ones(1,p.number_iterations)*p.gradthreshold;
            str_p = ['Applied per slice along axis 1, number of iterations: ' num2str(p.number_iterations,'%i') ', gradient thershold for each iteration: ' num2str(Gradient_threshold)];
            slice_ = imdiffusefilt(slice_,'Connectivity',p.connectivity,'ConductionMethod',p.conduction,'GradientThreshold',Gradient_threshold,'NumberOfIterations',p.number_iterations);
        end
        M(x,:,:) = slice_;        
    end

elseif strcmp(p.apply,'Slice per slice, along axis 2')
    for y=1:1:sz(2)
        slice_ = squeeze(M(:,y,:));
        if p.estimategradientiteration
            str_p = 'Applied per slice along axis 2, auto estimation of number of iterations and gradient threshold';
            [Gradient_threshold,number_iteration] = imdiffuseest(slice_,'Connectivity',p.connectivity,'ConductionMethod',p.conduction);
            slice_ = imdiffusefilt(slice_,'Connectivity',p.connectivity,'ConductionMethod',p.conduction,'GradientThreshold',Gradient_threshold,'NumberOfIterations',number_iteration);
        else
            Gradient_threshold = ones(1,p.number_iterations)*p.gradthreshold;
            str_p = ['Applied per slice along axis 2, number of iterations: ' num2str(p.number_iterations,'%i') ', gradient thershold for each iteration: ' num2str(Gradient_threshold)];
            slice_ = imdiffusefilt(slice_,'Connectivity',p.connectivity,'ConductionMethod',p.conduction,'GradientThreshold',Gradient_threshold,'NumberOfIterations',p.number_iterations);
        end
        M(:,y,:) = slice_;        
    end    

elseif strcmp(p.apply,'Slice per slice, along axis 3')
    for z=1:1:sz(3)
        if p.estimategradientiteration
            str_p = 'Applied per slice along axis 3, auto estimation of number of iterations and gradient threshold';
            [Gradient_threshold,number_iteration] = imdiffuseest(M(:,:,z),'Connectivity',p.connectivity,'ConductionMethod',p.conduction);
            M(:,:,z) = imdiffusefilt(M(:,:,z),'Connectivity',p.connectivity,'ConductionMethod',p.conduction,'GradientThreshold',Gradient_threshold,'NumberOfIterations',number_iteration);
        else
            Gradient_threshold = ones(1,p.number_iterations)*p.gradthreshold;
            str_p = ['Applied per slice along axis 3, number of iterations: ' num2str(p.number_iterations,'%i') ', gradient thershold for each iteration: ' num2str(Gradient_threshold)];
            M(:,:,z) = imdiffusefilt(M(:,:,z),'Connectivity',p.connectivity,'ConductionMethod',p.conduction,'GradientThreshold',Gradient_threshold,'NumberOfIterations',p.number_iterations);
        end        
    end
end

str_p =['Connectivity: ' p.connectivity ', Conduction method:' p.conduction ', ' str_p];

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