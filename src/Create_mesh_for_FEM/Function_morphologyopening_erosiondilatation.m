function [M] = Function_morphologyopening_erosiondilatation(M,background_code,erosion_distance,dilatation_distance,distance_method, iteration, approach, preserve_background)
% Morphology opening (erosion and dilation) to simplify binary image surface
% BW = Dilation( Erosion (BW) )

tolerance = 0.1;

for current_iteration = 1:1:iteration
    sz = size(M);
    Ma = zeros(sz)+background_code;
    phases = unique(M);
    phases(phases==background_code)=[];
    n_phase = length(phases);
    
    if strcmp(approach,'Solid phase per solid phase')
        for k=1:1:n_phase % Loop over all phase
            phase = phases(k);
            idx = find(M==phase);
            % Work on a subdomain to speed up calculation
            [I1,I2,I3] = ind2sub(sz,idx);
            x_min = min(I1); x_max = max(I1);
            y_min = min(I2); y_max = max(I2);
            z_min = min(I3); z_max = max(I3);
            subdomain = M(x_min:x_max, y_min:y_max, z_min:z_max);
            % Binarize the subdomain
            subdomain_BW = zeros(size(subdomain));
            idx = find(subdomain==phase);
            subdomain_BW(idx)=1;
            % Erosion dilatation applied on the subdomain to simplify interface
            [subdomain_BW] = Function_morphologyopening_erosiondilatation_algorithm(subdomain_BW,erosion_distance,tolerance,dilatation_distance,distance_method);
            % Re-insert in the full domain
            tmp = zeros(sz);
            tmp(x_min:x_max, y_min:y_max, z_min:z_max)=subdomain_BW;
            idx = find(tmp==1);
            Ma(idx)=phase;
        end
        
    elseif strcmp(approach,'All solid phases in one step')
        BW = zeros(sz);
        BW(M~=background_code)=1;
        [BW] = Function_morphologyopening_erosiondilatation_algorithm(BW,erosion_distance,tolerance,dilatation_distance,distance_method);
        cond_1 = BW == 1; % Within solid phase
        for k=1:1:n_phase % Loop over all phase
            phase = phases(k);
            cond_2 = M==phase;
            idx = find(cond_1+cond_2==2);
            Ma(idx) = phase;
        end
        cond_2 = Ma==background_code;
        idx = find(cond_1+cond_2==2);
        if ~isempty(idx)
            Edges_=zeros(sz);
            Edges_(idx) = 1;
            [~,idy] = bwdist(~Edges_,'Euclidean');
            for k=1:1:numel(idy)
                Ma(k)=M(idy(k));
            end
            Ma(BW==background_code)=background_code;
        end
        
    end
    
    if preserve_background % Assign background voxel to the nearest phase
        backgdound_location = find(M==background_code);
        BW = zeros(sz);
        BW(Ma~=background_code)=1;
        [~,idx] = bwdist(BW,'Euclidean');
        M=Ma;
        for k=1:1:numel(idx)
            M(k)=M(idx(k));
        end
        M(backgdound_location)=background_code;
    else
        M=Ma;
    end
    
end


end

