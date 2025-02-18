function [M,newtype,foo] = fct_calculatecluster(M,p)

newtype = 'Segmented (instance)';
foo = 1;

sz = size(M);
dimension = length(sz);

if ~p.background_is_edges
    idbackground = find(M==p.backgroundlabel);
end

BW = zeros(sz);
BW(M==p.label)=1;

if dimension == 2
    if strcmp(p.connectivity_rule_2D,'edges (4)')
        conn = 4;        
    elseif strcmp(p.connectivity_rule_2D,'edges or corners (8)')
        conn = 8;
    end
    M = bwlabel(BW,conn); % 4,8
else
    if strcmp(p.connectivity_rule_3D,'faces (6)')
        conn = 6;        
    elseif strcmp(p.connectivity_rule_3D,'faces or edges (18)')
        conn = 18;
    elseif strcmp(p.connectivity_rule_3D,'faces or edges or corners (26)')
        conn = 26;        
    end
    M = bwlabeln(BW,conn); % 6, 18, 26
end

if ~p.background_is_edges
    M = M+1;
    M(idbackground)=0;
end

[M] = fct_intconvert(M);

end