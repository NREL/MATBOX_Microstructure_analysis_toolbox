function [node, face, elem, subdomain] = function_iso2mesh_from_array(array_, options)

%% OPTIONS
method_surfacemesh = 'cgalmesh';
opt.radbound = options.radbound;
opt.distbound = options.distbound;
isovalues_ = unique(array_);

smooth_method = options.method_surfacemesh;
iteration_smoothing = options.iteration_smoothing;
useralpha = options.useralpha;
if strcmp('smooth_method','laplacianhc')
    userbeta = options.userbeta;
end

keepratio = options.keepratio;
maxvol = options.maxvol;

%% SURFACE MESH
tmp=uint8(array_); % Convert in unsigned integer 8bits
[nodes_coordinates,triangular_faces,regions,holes]=v2s(tmp,isovalues_,opt,method_surfacemesh);

%% SMOOTHING
conn=meshconn(triangular_faces(:,1:3),size(nodes_coordinates,1));
mask_=[];
smoothed_nodes_coordinates=smoothsurf(nodes_coordinates,mask_,conn,iteration_smoothing,useralpha,smooth_method);

%% REPAIRING MESH
do_repair=false;
if do_repair
    % This should prevent self-intersecting elements and fill holes
    %[smoothed_nodes_coordinates_r,triangular_faces_r]=meshcheckrepair(smoothed_nodes_coordinates,triangular_faces,'meshfix'); % 'meshfix' breaks the mesh
    [smoothed_nodes_coordinates_r,triangular_faces_r]=meshcheckrepair(smoothed_nodes_coordinates,triangular_faces,'deep'); % 
end

%% VOLUMETRIC MESH
Domain_size = size(array_);
x0=0; y0=0; z0=0;
x1=Domain_size(1); y1=Domain_size(2); z1=Domain_size(3);
if do_repair
    [node,elem,face]=surf2mesh(smoothed_nodes_coordinates_r,triangular_faces_r,[x0-1 y0-1 z0-1],[x1+1 y1+1 z1+1],keepratio,maxvol,regions);
else
    [node,elem,face]=surf2mesh(smoothed_nodes_coordinates,triangular_faces,[x0-1 y0-1 z0-1],[x1+1 y1+1 z1+1],keepratio,maxvol,regions);
end

%% SELECT RELEVANT COLUMNS
node = node(:,1:3);
subdomain = elem(:,5);
elem = elem(:,1:4);
face = face(:,1:3);

%% MESH CORRECTION
% Iso2mesh tends to lost the correct cell assignment for large volume
% This section needs more work
options.correction_misassigned_region = false;
options.correction_isolated_voxel = true;

number_region_expected = length(unique(array_));
allregion = unique(subdomain);
number_region_obtained = length(allregion);
disp(['Number of region expected: ' num2str(number_region_expected) ', number of region obtained: ' num2str(number_region_obtained)]);

if number_region_obtained>number_region_expected
    if options.correction_misassigned_region
        subdomain_id= unique(subdomain);
        tmp = zeros(size(subdomain));
        for k=1:1:length(subdomain_id)
            idx=find(subdomain==subdomain_id(k));
            % Cell centroid coordinate
            cell_centroid_x = ( node(elem(idx,1),1) + node(elem(idx,2),1) + node(elem(idx,3),1) + node(elem(idx,4),1) )/4;
            cell_centroid_y = ( node(elem(idx,1),2) + node(elem(idx,2),2) + node(elem(idx,3),2) + node(elem(idx,4),2) )/4;
            cell_centroid_z = ( node(elem(idx,1),3) + node(elem(idx,2),3) + node(elem(idx,3),3) + node(elem(idx,4),3) )/4;
            % Cell centroid array index
            cell_centroid_x = floor(cell_centroid_x)+1;
            cell_centroid_y = floor(cell_centroid_y)+1;
            cell_centroid_z = floor(cell_centroid_z)+1;
            % Evalue array
            ida = sub2ind(size(array_), cell_centroid_x, cell_centroid_y, cell_centroid_z);
            val = array_(ida);
            
            % Assign each cells of this subdomain to the array value it corresponds
            % tmp(idx)=val;
            
            % Assign all cells of this subdomain to the array value it corresponds the most
            u_val = unique(val);
            count_region = histc(val, u_val);
            arrayid = find(count_region==max(count_region));
            new_id = u_val(arrayid(1));
            tmp(idx)=new_id;
        end
        subdomain = tmp; clear tmp;
        number_region_obtained = length(unique(subdomain));
        disp(['Number of region expected: ' num2str(number_region_expected) ', number of region obtained after correction: ' num2str(number_region_obtained)]);
    end
    
    if options.correction_isolated_voxel
        count_region = histc(subdomain, allregion);
        id = find(count_region==min(count_region));
        incorrect_region = allregion(id(1));
        incorrect_region
        %incorrect_region = 0; % Unassigned cells are set to 0 by Iso2mesh by default
        tmp = unique(subdomain);
        tmp(tmp==incorrect_region)=[];
        other_region_max = max(tmp);
        other_region_min = min(tmp);
        for pass_=1:1:2
            incorrect_cell=1e9; % Initialize
            while incorrect_cell>0
                for other_region=other_region_max:-1:other_region_min
                    if other_region~=incorrect_region
                        change_=1e9; % Initialization
                        previous_incorrect_cell = incorrect_cell;
                        while change_~=0
                            
                            % All node that belong to the wrong zone
                            idx_wrongregion=find(subdomain==incorrect_region);
                            all_node_wrongregion=[elem(idx_wrongregion,1); elem(idx_wrongregion,2); elem(idx_wrongregion,3); elem(idx_wrongregion,4)];
                            all_node_wrongregion=unique(all_node_wrongregion);
                            
                            % Check if they belong to other zone
                            idx_otherregion=find(subdomain==other_region);
                            all_node_otherregion=[elem(idx_otherregion,1); elem(idx_otherregion,2); elem(idx_otherregion,3); elem(idx_otherregion,4)];
                            all_node_otherregion=unique(all_node_otherregion);
                            
                            % Node of wrong region that belong also to current other region
                            idx_tmp = ismember(all_node_wrongregion,all_node_otherregion);
                            idx_tmp = find(idx_tmp==1);
                            index_node_WrongInOther = all_node_wrongregion(idx_tmp);
                            
                            index_node_WrongInOther_a = ismember(elem(:,1),index_node_WrongInOther);
                            index_node_WrongInOther_b = ismember(elem(:,2),index_node_WrongInOther);
                            index_node_WrongInOther_c = ismember(elem(:,3),index_node_WrongInOther);
                            index_node_WrongInOther_d = ismember(elem(:,4),index_node_WrongInOther);
                            idx_a = find(index_node_WrongInOther_a==1);
                            idx_b = find(index_node_WrongInOther_b==1);
                            idx_c = find(index_node_WrongInOther_c==1);
                            idx_d = find(index_node_WrongInOther_d==1);
                            idx_all=[idx_a; idx_b; idx_c; idx_d];
                            idx_all=unique(idx_all);
                            
                            tmp = subdomain(idx_all);
                            tmp2 = find(tmp(:,1)==incorrect_region);
                            idx_all_ = idx_all(tmp2);
                            
                            incorrect_cell=length(idx_all_);
                            fprintf ('- Pass %i. Check incorrect region with region %i, number of incorrect cell %i\n',pass_,other_region,incorrect_cell);
                            
                            change_ = incorrect_cell - previous_incorrect_cell;
                            previous_incorrect_cell = incorrect_cell;
                            
                            subdomain(idx_all_)=other_region;
                        end
                    end
                end
            end
        end
        number_region_obtained = length(unique(subdomain));
        disp(['Number of region expected: ' num2str(number_region_expected) ', number of region obtained after correction: ' num2str(number_region_obtained)]);
    end
    
end

%% FIX GROUP CONNECTIVITY
% Group should be contiguous. Check connectivty and reassign disconnected clusters


%% FIX VERTICE-VERTICE AND EDGE-EDGE ONLY CONNECTIONS
% Cell in contact with the rest of the subdomain with only one vertice or only one edge (line) will be reassigned to the subdomain for which they share most of their faces



end