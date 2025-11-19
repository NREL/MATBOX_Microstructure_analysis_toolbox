function [M,newtype,foo] = fct_globalerosiondilation(M,p)

sz = size(M);
foo=[];
tolerance = 0.01;

dimension = length(sz);

if p.zeroes_on_edges
    d = max([p.erosiondistance, p.dilationdistance])+2;
    sz = sz+2*d;
    Moriginal = M;
    M = zeros(sz);
    % if dimension==2
    %     M(2:end-1,2:end-1) = Moriginal;
    % else
    %     M(2:end-1,2:end-1,2:end-1) = Moriginal;
    % end

    if dimension==2
        M(d+1:end-d,d+1:end-d) = Moriginal;
    else
        M(d+1:end-d,d+1:end-d,d+1:end-d) = Moriginal;
    end

    clear Moriginal;
end

if strcmp(p.applyon,'A single label')
    newtype = 'same';
    BW = zeros(sz,'uint8');
    idx = M==p.labeltoapply;
    BW(idx)=1;

    for iter=1:1:p.num_iteration
        if strcmp(p.order,'Opening (erosion then dilation)')
            if p.erosiondistance>0
                distance_map = bwdist(~BW,p.distancedefinition);
                BW(distance_map<=p.erosiondistance+tolerance)=0;
            end
            if p.dilationdistance>0
                distance_map = bwdist(BW,p.distancedefinition);
                BW(distance_map<=p.dilationdistance+tolerance)=1;
            end
        elseif strcmp(p.order,'Closing (dilation then erosion)')
            if p.dilationdistance>0
                distance_map = bwdist(BW,p.distancedefinition);
                BW(distance_map<=p.dilationdistance+tolerance)=1;
            end
            if p.erosiondistance>0
                distance_map = bwdist(~BW,p.distancedefinition);
                BW(distance_map<=p.erosiondistance+tolerance)=0;
            end
        end
    end

    M(idx)=p.background; % Remove
    n_nonbackgroundlabel = length(unique(M))-1;
    if p.dilationcanoverwritte || n_nonbackgroundlabel==0 
        M(BW==1)=p.labeltoapply; % Add and overwritte other non-background label(s)
    else
        idxother = M~=p.background;
        Msav = M;
        M(BW==1)=p.labeltoapply; % Add...
        M(idxother)=Msav(idxother); % ... but do not overwritte other non-background label(s)
    end

else

    newtype = 'Segmented (phase)';
    % Particle per particle to avoid merging
    unique_labels = unique(M);
    unique_labels(unique_labels==p.background)=[];
    n_labels = length(unique_labels);
    BW = zeros(sz);
    for k=1:1:n_labels
        BWtmp = zeros(sz);
        BWtmp(M==unique_labels(k))=1;

        for iter=1:1:p.num_iteration
            if strcmp(p.order,'Opening (erosion then dilation)')
                if p.erosiondistance>0
                    distance_map = bwdist(~BWtmp,p.distancedefinition);
                    BWtmp(distance_map<=p.erosiondistance+tolerance)=0;
                end
                if p.dilationdistance>0
                    distance_map = bwdist(BWtmp,p.distancedefinition);
                    BWtmp(distance_map<=p.dilationdistance+tolerance)=1;
                end
            elseif strcmp(p.order,'Closing (dilation then erosion)')
                if p.dilationdistance>0
                    distance_map = bwdist(BWtmp,p.distancedefinition);
                    BWtmp(distance_map<=p.dilationdistance+tolerance)=1;
                end
                if p.erosiondistance>0
                    distance_map = bwdist(~BWtmp,p.distancedefinition);
                    BWtmp(distance_map<=p.erosiondistance+tolerance)=0;
                end
            end
        end

        BW = BW + BWtmp;
    end
    BW(BW~=0)=1;
    M=BW;

end

if p.zeroes_on_edges
    if dimension==2
        M=M(d+1:end-d,d+1:end-d);
    else
        M=M(d+1:end-d,d+1:end-d,d+1:end-d);
    end    
    
    % if dimension==2
    %     M=M(2:end-1,2:end-1);
    % else
    %     M=M(2:end-1,2:end-1,2:end-1);
    % end
end

[M] = fct_intconvert(M);

end