function [node_,face_, elem_, subdomain_] = function_regularmesh_from_array(M,ids,cellchoice,do_face)

check_facedoublons = do_face; % % false to save RAM as optional step

number_phase = length(ids);
for kphase = 1:1:number_phase % Loop over phase
    BW = zeros(size(M)); % Initialize
    BW(M==ids(kphase))=1;
    [tmp.node_,tmp.facets, tmp.cells] = function_create_vertices_cell_from_array(BW,cellchoice,do_face);
    [number_cell, ~] = size(tmp.cells);
    tmp.subdomain_ = zeros(number_cell,1) + ids(kphase);
    % Concatenate
    if kphase==1
        vertices = tmp.node_(:,1:3);
        facets = tmp.facets;
        cells = tmp.cells;
        subdomain_ = tmp.subdomain_;
    else
        [indexstart,~] = size(vertices);
        vertices = [vertices; tmp.node_(:,1:3)];
        facets = [facets; tmp.facets+indexstart];
        cells = [cells; tmp.cells+indexstart];
        subdomain_ = [subdomain_; tmp.subdomain_];
    end
end

if number_phase==1
    node_ = vertices;
    face_ = facets;
    elem_ = cells;
else
    % Remove vertices doublons using a pairing function
    % *10 so that x.5 floating coordinates are converted to x5 integers
    cantor_value = pairing_cantor_function(vertices(:,1)*10,vertices(:,2)*10);
    cantor_value = pairing_cantor_function(cantor_value,vertices(:,3)*10);
    [unique_cantor,ia,ic] = unique(cantor_value);
    
    % Update as consequence indexes in cells
    [number_cell, vertice_per_cell] = size(cells);
    elem_=zeros(number_cell, vertice_per_cell);
    node_ = vertices(ia,:);
    for line=1:1:number_cell
        for column=1:1:vertice_per_cell
            old = cells(line,column);
            elem_(line,column) = ic(old);
        end
    end
    
    % Update as consequence indexes in facets
    if check_facedoublons
        [number_facets, vertice_per_facets] = size(facets);
        face_=zeros(number_facets, vertice_per_facets);
        for line=1:1:number_facets
            for column=1:1:vertice_per_facets
                old = facets(line,column);
                face_(line,column) = ic(old);
            end
        end
        clear old % Cleaning memory
        % For multi phase, some facets are repeated twice and must be removed
        cantor_value = pairing_cantor_function(face_(:,1),face_(:,2));
        for kpairing=1:1:vertice_per_facets-2
            cantor_value = pairing_cantor_function(cantor_value,face_(:,kpairing+2));
        end
        [unique_cantor,ia,ic] = unique(cantor_value);
        face_ = face_(ia,:); % BUG
    else
        face_ = facets;
    end
end

    function k3 = pairing_cantor_function(k1,k2)
        k3 = 0.5.*(k1+k2).*(k1+k2+1)+k2;
    end

end