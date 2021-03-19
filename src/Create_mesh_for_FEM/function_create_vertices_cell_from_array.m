function [node_,face_, elem_] = function_create_vertices_cell_from_array(BW,cellchoice,do_face)
% [node_,elem_] = function_create_vertices_cell_from_array(BW,cellchoice,indexstart)
% Inputs:
% - BW a 3D binary array with BW(i,j,k)=1 if voxel belong to the investigated phase
% - cellchoice: '5_tetrahedron_per_voxel' or '24_tetrahedron_per_voxel'
% - do face: true or false. Facets are not required for re-generating the
% mesh, but are required for visualization
% Outputs:
% - node_: n*3 array: x y z coordinates of vertices
% - face_: m*3 array: a b c index of vertices of triangle face
% - elem_: k*4 array: a b c d index of vertices of tetrahedron cells

number_voxelphase = sum(sum(sum(BW==1)));
if strcmp(cellchoice,'5 tetrahedrons')
    cell_per_voxel = 5;
    vertice_per_cell = 4;
elseif strcmp(cellchoice,'6 tetrahedrons')
    cell_per_voxel = 6;
    vertice_per_cell = 4;
elseif strcmp(cellchoice,'24 tetrahedrons')
    cell_per_voxel = 24;
    vertice_per_cell = 4;
elseif strcmp(cellchoice,'1 hexahedron')
    cell_per_voxel = 1;
    vertice_per_cell = 8;
end
elem_ = zeros(cell_per_voxel*number_voxelphase,vertice_per_cell);

% Initialize
[number_cell, ~] = size(elem_);
elem_=zeros(number_cell, vertice_per_cell);

% Find voxels
domain_size = size(BW);
idx = find(BW==1);
[IX,IY,IZ] = ind2sub(domain_size,idx);

% Voxels min-max coordinate
x_min = IX-1; x_max = IX; x_mean = (x_min+x_max)/2;
y_min = IY-1; y_max = IY; y_mean = (y_min+y_max)/2;
z_min = IZ-1; z_max = IZ; z_mean = (z_min+z_max)/2;

% x_min = IX-1+0.5; x_max = IX+0.5; x_mean = (x_min+x_max)/2;
% y_min = IY-1+0.5; y_max = IY+0.5; y_mean = (y_min+y_max)/2;
% z_min = IZ-1+0.5; z_max = IZ+0.5; z_mean = (z_min+z_max)/2;

if strcmp(cellchoice,'5 tetrahedrons')
    % Vertices
    vertex_A = [x_min y_min z_min (1:1:number_voxelphase)'];
    vertex_B = [x_min y_max z_min (1:1:number_voxelphase)'+number_voxelphase];
    vertex_C = [x_max y_max z_min (1:1:number_voxelphase)'+2*number_voxelphase];
    vertex_D = [x_max y_min z_min (1:1:number_voxelphase)'+3*number_voxelphase];
    vertex_E = [x_min y_min z_max (1:1:number_voxelphase)'+4*number_voxelphase];
    vertex_F = [x_min y_max z_max (1:1:number_voxelphase)'+5*number_voxelphase];
    vertex_G = [x_max y_max z_max (1:1:number_voxelphase)'+6*number_voxelphase];
    vertex_H = [x_max y_min z_max (1:1:number_voxelphase)'+7*number_voxelphase];
    vertices = [vertex_A; vertex_B; vertex_C; vertex_D; vertex_E; vertex_F; vertex_G; vertex_H];
    % Cells
    tet_1 = [vertex_A(:,4) vertex_B(:,4) vertex_C(:,4) vertex_F(:,4)];
    tet_2 = [vertex_A(:,4) vertex_D(:,4) vertex_H(:,4) vertex_C(:,4)];
    tet_3 = [vertex_C(:,4) vertex_F(:,4) vertex_G(:,4) vertex_H(:,4)];
    tet_4 = [vertex_A(:,4) vertex_F(:,4) vertex_H(:,4) vertex_C(:,4)]; % Regular central tetrahedron
    tet_5 = [vertex_A(:,4) vertex_F(:,4) vertex_H(:,4) vertex_E(:,4)];
    if ~do_face % Memory managment
        clear vertex_A vertex_B vertex_C vertex_D vertex_E vertex_F vertex_G vertex_H
    end    
    cells = [tet_1; tet_2; tet_3; tet_4; tet_5];
    clear tet_1 tet_2 tet_3 tet_4 tet_5 % Memory managment
    % Faces
    if do_face
        face_ABF = [vertex_A(:,4) vertex_B(:,4) vertex_F(:,4)];
        face_BFC = [vertex_B(:,4) vertex_F(:,4) vertex_C(:,4)];
        face_FCA = [vertex_F(:,4) vertex_C(:,4) vertex_A(:,4)];
        face_BCA = [vertex_B(:,4) vertex_C(:,4) vertex_A(:,4)];
        
        face_AFH = [vertex_A(:,4) vertex_F(:,4) vertex_H(:,4)];
        face_AFE = [vertex_A(:,4) vertex_F(:,4) vertex_E(:,4)];
        face_FHE = [vertex_F(:,4) vertex_H(:,4) vertex_E(:,4)];
        face_AEH = [vertex_A(:,4) vertex_E(:,4) vertex_H(:,4)];
        
        face_FCG = [vertex_F(:,4) vertex_C(:,4) vertex_G(:,4)];
        face_GHC = [vertex_G(:,4) vertex_H(:,4) vertex_C(:,4)];
        face_FCH = [vertex_F(:,4) vertex_C(:,4) vertex_H(:,4)];
        face_FHG = [vertex_F(:,4) vertex_H(:,4) vertex_G(:,4)];
        
        face_ADC = [vertex_A(:,4) vertex_D(:,4) vertex_C(:,4)];
        face_CHD = [vertex_C(:,4) vertex_H(:,4) vertex_D(:,4)];
        face_ADH = [vertex_A(:,4) vertex_D(:,4) vertex_H(:,4)];
        face_AHC = [vertex_A(:,4) vertex_H(:,4) vertex_C(:,4)];
        clear vertex_A vertex_B vertex_C vertex_D vertex_E vertex_F vertex_G vertex_H
        facets = [face_ABF; face_BFC; face_FCA; face_BCA; face_AFH; face_AFE; face_FHE; face_AEH; face_FCG; face_GHC; face_FCH; face_FHG; face_ADC; face_CHD; face_ADH; face_AHC];
        clear face_ABF face_BFC face_FCA face_BCA face_AFH face_AFE face_FHE face_AEH face_FCG face_GHC face_FCH face_FHG face_ADC face_CHD face_ADH face_AHC
    else
        facets=[];
    end
    
elseif strcmp(cellchoice,'6 tetrahedrons')
    % Vertices
    vertex_A = [x_min y_min z_min (1:1:number_voxelphase)'];
    vertex_B = [x_min y_max z_min (1:1:number_voxelphase)'+number_voxelphase];
    vertex_C = [x_max y_max z_min (1:1:number_voxelphase)'+2*number_voxelphase];
    vertex_D = [x_max y_min z_min (1:1:number_voxelphase)'+3*number_voxelphase];
    vertex_E = [x_min y_min z_max (1:1:number_voxelphase)'+4*number_voxelphase];
    vertex_F = [x_min y_max z_max (1:1:number_voxelphase)'+5*number_voxelphase];
    vertex_G = [x_max y_max z_max (1:1:number_voxelphase)'+6*number_voxelphase];
    vertex_H = [x_max y_min z_max (1:1:number_voxelphase)'+7*number_voxelphase];
    vertices = [vertex_A; vertex_B; vertex_C; vertex_D; vertex_E; vertex_F; vertex_G; vertex_H];
    
    % Cells
    tet_1 = [vertex_A(:,4) vertex_B(:,4) vertex_F(:,4) vertex_G(:,4)];
    tet_2 = [vertex_A(:,4) vertex_B(:,4) vertex_C(:,4) vertex_G(:,4)];
    tet_3 = [vertex_A(:,4) vertex_F(:,4) vertex_G(:,4) vertex_E(:,4)];
    tet_4 = [vertex_A(:,4) vertex_C(:,4) vertex_G(:,4) vertex_D(:,4)];
    tet_5 = [vertex_A(:,4) vertex_E(:,4) vertex_G(:,4) vertex_H(:,4)];
    tet_6 = [vertex_A(:,4) vertex_D(:,4) vertex_G(:,4) vertex_H(:,4)];
    if ~do_face % Memory managment
        clear vertex_A vertex_B vertex_C vertex_D vertex_E vertex_F vertex_G vertex_H
    end
    cells = [tet_1; tet_2; tet_3; tet_4; tet_5; tet_6];
    clear tet_1 tet_2 tet_3 tet_4 tet_5 tet_6
        
    % Faces
    if do_face
        face_ABG = [vertex_A(:,4) vertex_B(:,4) vertex_G(:,4)];
        face_ABF = [vertex_A(:,4) vertex_B(:,4) vertex_F(:,4)];
        face_BFG = [vertex_B(:,4) vertex_F(:,4) vertex_G(:,4)];
        face_AFG = [vertex_A(:,4) vertex_F(:,4) vertex_G(:,4)];
        
        face_ABC = [vertex_A(:,4) vertex_B(:,4) vertex_C(:,4)];
        face_BCG = [vertex_B(:,4) vertex_C(:,4) vertex_G(:,4)];
        face_ACG = [vertex_A(:,4) vertex_C(:,4) vertex_G(:,4)];
        
        face_AEF = [vertex_A(:,4) vertex_E(:,4) vertex_F(:,4)];
        face_EFG = [vertex_E(:,4) vertex_F(:,4) vertex_G(:,4)];
        face_AEG = [vertex_A(:,4) vertex_E(:,4) vertex_G(:,4)];
        
        face_ADC = [vertex_A(:,4) vertex_D(:,4) vertex_C(:,4)];
        face_DCG = [vertex_D(:,4) vertex_C(:,4) vertex_G(:,4)];
        face_ADG = [vertex_A(:,4) vertex_D(:,4) vertex_G(:,4)];
        
        face_AEH = [vertex_A(:,4) vertex_E(:,4) vertex_H(:,4)];
        face_AHG = [vertex_A(:,4) vertex_H(:,4) vertex_G(:,4)];
        face_EHG = [vertex_E(:,4) vertex_H(:,4) vertex_G(:,4)];
        
        face_AHD = [vertex_A(:,4) vertex_H(:,4) vertex_D(:,4)];
        face_HGD = [vertex_H(:,4) vertex_G(:,4) vertex_D(:,4)];
        
        clear vertex_A vertex_B vertex_C vertex_D vertex_E vertex_F vertex_G vertex_H
        facets = [face_ABG; face_ABF; face_BFG; face_AFG; face_ABC; face_BCG; face_ACG; face_AEF; face_EFG; face_AEG; face_ADC; face_DCG; face_ADG; face_AEH; face_AHG; face_EHG; face_AHD; face_HGD];
        clear face_ABG face_ABF face_BFG face_AFG face_ABC face_BCG face_ACG face_AEF face_EFG face_AEG face_ADC face_DCG face_ADG face_AEH face_AHG face_EHG face_AHD face_HGD
    else 
        facets=[];
    end
    
elseif strcmp(cellchoice,'24 tetrahedrons')
    % Vertices
    % Voxel corners
    vertex_A = [x_min y_min z_min (1:1:number_voxelphase)'];
    vertex_B = [x_min y_max z_min (1:1:number_voxelphase)'+number_voxelphase];
    vertex_C = [x_max y_max z_min (1:1:number_voxelphase)'+2*number_voxelphase];
    vertex_D = [x_max y_min z_min (1:1:number_voxelphase)'+3*number_voxelphase];
    vertex_E = [x_min y_min z_max (1:1:number_voxelphase)'+4*number_voxelphase];
    vertex_F = [x_min y_max z_max (1:1:number_voxelphase)'+5*number_voxelphase];
    vertex_G = [x_max y_max z_max (1:1:number_voxelphase)'+6*number_voxelphase];
    vertex_H = [x_max y_min z_max (1:1:number_voxelphase)'+7*number_voxelphase];
    % Face center
    vertex_ABCD = [x_mean y_mean z_min (1:1:number_voxelphase)'+8*number_voxelphase];
    vertex_EFGH = [x_mean y_mean z_max (1:1:number_voxelphase)'+9*number_voxelphase];
    vertex_ADHE = [x_mean y_min z_mean (1:1:number_voxelphase)'+10*number_voxelphase];
    vertex_BCGF = [x_mean y_max z_mean (1:1:number_voxelphase)'+11*number_voxelphase];
    vertex_AEFB = [x_min y_mean z_mean (1:1:number_voxelphase)'+12*number_voxelphase];
    vertex_DCGH = [x_max y_mean z_mean (1:1:number_voxelphase)'+13*number_voxelphase];
    % Voxel centroid
    vertex_O = [x_mean y_mean z_mean (1:1:number_voxelphase)'+14*number_voxelphase];
    vertices = [vertex_A; vertex_B; vertex_C; vertex_D; vertex_E; vertex_F; vertex_G; vertex_H; vertex_ABCD; vertex_EFGH; vertex_ADHE; vertex_BCGF; vertex_AEFB; vertex_DCGH; vertex_O];
    
    % Cells 
    tet_1 = [vertex_A(:,4) vertex_B(:,4) vertex_ABCD(:,4) vertex_O(:,4)];
    tet_2 = [vertex_B(:,4) vertex_C(:,4) vertex_ABCD(:,4) vertex_O(:,4)];
    tet_3 = [vertex_C(:,4) vertex_D(:,4) vertex_ABCD(:,4) vertex_O(:,4)];
    tet_4 = [vertex_D(:,4) vertex_A(:,4) vertex_ABCD(:,4) vertex_O(:,4)];
    tet_5 = [vertex_E(:,4) vertex_F(:,4) vertex_EFGH(:,4) vertex_O(:,4)];
    tet_6 = [vertex_F(:,4) vertex_G(:,4) vertex_EFGH(:,4) vertex_O(:,4)];
    tet_7 = [vertex_G(:,4) vertex_H(:,4) vertex_EFGH(:,4) vertex_O(:,4)];
    tet_8 = [vertex_H(:,4) vertex_E(:,4) vertex_EFGH(:,4) vertex_O(:,4)];
    tet_9  = [vertex_A(:,4) vertex_D(:,4) vertex_ADHE(:,4) vertex_O(:,4)];
    tet_10 = [vertex_D(:,4) vertex_H(:,4) vertex_ADHE(:,4) vertex_O(:,4)];
    tet_11 = [vertex_H(:,4) vertex_E(:,4) vertex_ADHE(:,4) vertex_O(:,4)];
    tet_12 = [vertex_E(:,4) vertex_A(:,4) vertex_ADHE(:,4) vertex_O(:,4)];
    tet_13  = [vertex_B(:,4) vertex_C(:,4) vertex_BCGF(:,4) vertex_O(:,4)];
    tet_14 = [vertex_C(:,4) vertex_G(:,4) vertex_BCGF(:,4) vertex_O(:,4)];
    tet_15 = [vertex_G(:,4) vertex_F(:,4) vertex_BCGF(:,4) vertex_O(:,4)];
    tet_16 = [vertex_F(:,4) vertex_B(:,4) vertex_BCGF(:,4) vertex_O(:,4)];
    tet_17  = [vertex_A(:,4) vertex_E(:,4) vertex_AEFB(:,4) vertex_O(:,4)];
    tet_18 = [vertex_E(:,4) vertex_F(:,4) vertex_AEFB(:,4) vertex_O(:,4)];
    tet_19 = [vertex_F(:,4) vertex_B(:,4) vertex_AEFB(:,4) vertex_O(:,4)];
    tet_20 = [vertex_B(:,4) vertex_A(:,4) vertex_AEFB(:,4) vertex_O(:,4)];
    tet_21  = [vertex_D(:,4) vertex_C(:,4) vertex_DCGH(:,4) vertex_O(:,4)];
    tet_22 = [vertex_C(:,4) vertex_G(:,4) vertex_DCGH(:,4) vertex_O(:,4)];
    tet_23 = [vertex_G(:,4) vertex_H(:,4) vertex_DCGH(:,4) vertex_O(:,4)];
    tet_24 = [vertex_H(:,4) vertex_D(:,4) vertex_DCGH(:,4) vertex_O(:,4)];
    if ~do_face % Memory managment
        clear vertex_A vertex_B vertex_C vertex_D vertex_E vertex_F vertex_G vertex_H vertex_ABCD vertex_EFGH vertex_ADHE vertex_BCGF vertex_AEFB vertex_DCGH vertex_O
    end    
    cells = [tet_1; tet_2; tet_3; tet_4; tet_5; tet_6; tet_7; tet_8; tet_9; tet_10; tet_11; tet_12; tet_13; tet_14; tet_15; tet_16; tet_17; tet_18; tet_19; tet_20; tet_21; tet_22; tet_23; tet_24];
    clear tet_1 tet_2 tet_3 tet_4 tet_5 tet_6 tet_7 tet_8 tet_9 tet_10 tet_11 tet_12 tet_13 tet_14 tet_15 tet_16 tet_17 tet_18 tet_19 tet_20 tet_21 tet_22 tet_23 tet_24
    
    % Faces
    if do_face
        face_ABO = [vertex_A(:,4) vertex_B(:,4) vertex_O(:,4)];
        face_BC1O = [vertex_B(:,4) vertex_AEFB(:,4) vertex_O(:,4)];
        face_AC1O = [vertex_A(:,4) vertex_AEFB(:,4) vertex_O(:,4)];
        face_ABC1 = [vertex_A(:,4) vertex_B(:,4) vertex_AEFB(:,4)];
        
        face_BFO = [vertex_B(:,4) vertex_F(:,4) vertex_O(:,4)];
        face_FC1O = [vertex_F(:,4) vertex_AEFB(:,4) vertex_O(:,4)];
        face_BFC1 = [vertex_B(:,4) vertex_F(:,4) vertex_AEFB(:,4)];
        
        face_FEO = [vertex_F(:,4) vertex_E(:,4) vertex_O(:,4)];
        face_EC1O = [vertex_E(:,4) vertex_AEFB(:,4) vertex_O(:,4)];
        face_FEC1 = [vertex_F(:,4) vertex_E(:,4) vertex_AEFB(:,4)];
        
        face_EAO = [vertex_E(:,4) vertex_A(:,4) vertex_O(:,4)];
        face_EAC1 = [vertex_E(:,4) vertex_A(:,4) vertex_AEFB(:,4)];
        
        
        face_BC2O = [vertex_B(:,4) vertex_BCGF(:,4) vertex_O(:,4)];
        face_FC2O = [vertex_F(:,4) vertex_BCGF(:,4) vertex_O(:,4)];
        face_BFC2 = [vertex_B(:,4) vertex_F(:,4) vertex_BCGF(:,4)];
        
        face_FGO = [vertex_F(:,4) vertex_G(:,4) vertex_O(:,4)];
        face_GC2O = [vertex_G(:,4) vertex_BCGF(:,4) vertex_O(:,4)];
        face_FGC2 = [vertex_F(:,4) vertex_G(:,4) vertex_BCGF(:,4)];
        
        face_GCO = [vertex_G(:,4) vertex_C(:,4) vertex_O(:,4)];
        face_CC2O = [vertex_C(:,4) vertex_BCGF(:,4) vertex_O(:,4)];
        face_GCC2 = [vertex_G(:,4) vertex_C(:,4) vertex_BCGF(:,4)];
        
        face_CBO = [vertex_C(:,4) vertex_B(:,4) vertex_O(:,4)];
        face_CBC2 = [vertex_C(:,4) vertex_B(:,4) vertex_BCGF(:,4)];
        
        
        face_FC3O = [vertex_F(:,4) vertex_EFGH(:,4) vertex_O(:,4)];
        face_GC3O = [vertex_G(:,4) vertex_EFGH(:,4) vertex_O(:,4)];
        face_FGC3 = [vertex_F(:,4) vertex_G(:,4) vertex_EFGH(:,4)];
        
        face_GHO = [vertex_G(:,4) vertex_H(:,4) vertex_O(:,4)];
        face_HC3O = [vertex_H(:,4) vertex_EFGH(:,4) vertex_O(:,4)];
        face_GHC3 = [vertex_G(:,4) vertex_H(:,4) vertex_EFGH(:,4)];
        
        face_HEO = [vertex_H(:,4) vertex_E(:,4) vertex_O(:,4)];
        face_EC3O = [vertex_E(:,4) vertex_EFGH(:,4) vertex_O(:,4)];
        face_HEC3 = [vertex_H(:,4) vertex_E(:,4) vertex_EFGH(:,4)];
        
        face_EFO = [vertex_E(:,4) vertex_F(:,4) vertex_O(:,4)];
        face_EFC3 = [vertex_E(:,4) vertex_F(:,4) vertex_EFGH(:,4)];
        
        
        face_CC4O = [vertex_C(:,4) vertex_DCGH(:,4) vertex_O(:,4)];
        face_GC4O = [vertex_G(:,4) vertex_DCGH(:,4) vertex_O(:,4)];
        face_CGC4 = [vertex_C(:,4) vertex_G(:,4) vertex_DCGH(:,4)];
        
        face_HC4O = [vertex_H(:,4) vertex_DCGH(:,4) vertex_O(:,4)];
        face_GHC4 = [vertex_G(:,4) vertex_H(:,4) vertex_DCGH(:,4)];
        
        face_HDO = [vertex_H(:,4) vertex_D(:,4) vertex_O(:,4)];
        face_DC4O = [vertex_D(:,4) vertex_DCGH(:,4) vertex_O(:,4)];
        face_HDC4 = [vertex_H(:,4) vertex_D(:,4) vertex_DCGH(:,4)];
        
        face_DCO = [vertex_D(:,4) vertex_C(:,4) vertex_O(:,4)];
        face_DCC4 = [vertex_D(:,4) vertex_C(:,4) vertex_DCGH(:,4)];
        
        
        face_AC5O = [vertex_A(:,4) vertex_ABCD(:,4) vertex_O(:,4)];
        face_BC5O = [vertex_B(:,4) vertex_ABCD(:,4) vertex_O(:,4)];
        face_ABC5 = [vertex_A(:,4) vertex_B(:,4) vertex_ABCD(:,4)];
        
        face_CC5O = [vertex_C(:,4) vertex_ABCD(:,4) vertex_O(:,4)];
        face_BCC5 = [vertex_B(:,4) vertex_C(:,4) vertex_ABCD(:,4)];
        
        face_DC5O = [vertex_D(:,4) vertex_ABCD(:,4) vertex_O(:,4)];
        face_CDC5 = [vertex_C(:,4) vertex_D(:,4) vertex_ABCD(:,4)];
        
        face_DAO = [vertex_D(:,4) vertex_A(:,4) vertex_O(:,4)];
        face_DAC5 = [vertex_D(:,4) vertex_A(:,4) vertex_ABCD(:,4)];
        
        
        face_AC6O = [vertex_A(:,4) vertex_ADHE(:,4) vertex_O(:,4)];
        face_EC6O = [vertex_E(:,4) vertex_ADHE(:,4) vertex_O(:,4)];
        face_AEC6 = [vertex_A(:,4) vertex_E(:,4) vertex_ADHE(:,4)];
        
        face_HC6O = [vertex_H(:,4) vertex_ADHE(:,4) vertex_O(:,4)];
        face_EHC6 = [vertex_E(:,4) vertex_H(:,4) vertex_ADHE(:,4)];
        
        face_DC6O = [vertex_D(:,4) vertex_ADHE(:,4) vertex_O(:,4)];
        face_HDC6 = [vertex_H(:,4) vertex_D(:,4) vertex_ADHE(:,4)];
        
        face_ADC6 = [vertex_A(:,4) vertex_D(:,4) vertex_ADHE(:,4)];
        
        clear vertex_A vertex_B vertex_C vertex_D vertex_E vertex_F vertex_G vertex_H vertex_ABCD vertex_EFGH vertex_ADHE vertex_BCGF vertex_AEFB vertex_DCGH vertex_O
        
        facets = [face_ABO; face_BC1O; face_AC1O; face_ABC1; face_BFO; face_FC1O; face_BFC1; face_FEO; face_EC1O; face_FEC1; face_EAO; face_EAC1];
        clear face_ABO face_BC1O face_AC1O face_ABC1 face_BFO face_FC1O face_BFC1 face_FEO face_EC1O face_FEC1 face_EAO face_EAC1
        facets = [facets; face_BC2O; face_FC2O; face_BFC2; face_FGO; face_GC2O; face_FGC2; face_GCO; face_CC2O; face_GCC2; face_CBO; face_CBC2];
        clear face_BC2O face_FC2O face_BFC2 face_FGO face_GC2O face_FGC2 face_GCO face_CC2O face_GCC2 face_CBO face_CBC
        facets = [facets; face_FC3O; face_GC3O; face_FGC3; face_GHO; face_HC3O; face_GHC3; face_HEO; face_EC3O; face_HEC3; face_EFO; face_EFC3];
        clear face_FC3O face_GC3O face_FGC3 face_GHO face_HC3O face_GHC3 face_HEO face_EC3O face_HEC3 face_EFO face_EFC3
        facets = [facets; face_CC4O; face_GC4O; face_CGC4; face_HC4O; face_GHC4; face_HDO; face_DC4O; face_HDC4; face_DCO; face_DCC4];
        clear face_CC4O face_GC4O face_CGC4 face_HC4O face_GHC4 face_HDO face_DC4O face_HDC4 face_DCO face_DCC4
        facets = [facets; face_AC5O; face_BC5O; face_ABC5; face_CC5O; face_BCC5; face_DC5O; face_CDC5; face_DAO; face_DAC5];
        clear face_AC5O face_BC5O face_ABC5 face_CC5O face_BCC5 face_DC5O face_CDC5 face_DAO face_DAC
        facets = [facets; face_AC6O; face_EC6O; face_AEC6; face_HC6O; face_EHC6; face_DC6O; face_HDC6; face_ADC6];
        clear face_AC6O face_EC6O face_AEC6 face_HC6O face_EHC6 face_DC6O face_HDC6 face_ADC6
    else
        facets=[];
    end
    
elseif strcmp(cellchoice,'1 hexahedron')
    % Vertices
    vertex_A = [x_min y_min z_min (1:1:number_voxelphase)'];
    vertex_B = [x_min y_max z_min (1:1:number_voxelphase)'+number_voxelphase];
    vertex_C = [x_max y_max z_min (1:1:number_voxelphase)'+2*number_voxelphase];
    vertex_D = [x_max y_min z_min (1:1:number_voxelphase)'+3*number_voxelphase];
    vertex_E = [x_min y_min z_max (1:1:number_voxelphase)'+4*number_voxelphase];
    vertex_F = [x_min y_max z_max (1:1:number_voxelphase)'+5*number_voxelphase];
    vertex_G = [x_max y_max z_max (1:1:number_voxelphase)'+6*number_voxelphase];
    vertex_H = [x_max y_min z_max (1:1:number_voxelphase)'+7*number_voxelphase];
    vertices = [vertex_A; vertex_B; vertex_C; vertex_D; vertex_E; vertex_F; vertex_G; vertex_H];
    % Cells
    hexa_1 = [vertex_A(:,4) vertex_B(:,4) vertex_C(:,4) vertex_D(:,4) vertex_E(:,4) vertex_F(:,4) vertex_G(:,4) vertex_H(:,4)];
    cells = [hexa_1];
    facets=[];
end


% Remove vertices doublons using a pairing function
% *10 so that x.5 floating coordinates are converted to x5 integers
cantor_value = pairing_cantor_function(vertices(:,1)*10,vertices(:,2)*10);
cantor_value = pairing_cantor_function(cantor_value,vertices(:,3)*10);
[unique_cantor,ia,ic] = unique(cantor_value);
node_ = vertices(ia,:);

% Update as consequence indexes in cells
for line=1:1:number_cell
    for column=1:1:vertice_per_cell
        old = cells(line,column);
        elem_(line,column) = ic(old);
    end
end

% Update as consequence indexes in facets
if do_face
    [number_facets, vertice_per_facets] = size(facets);
    face_=zeros(number_facets, vertice_per_facets);
    for line=1:1:number_facets
        for column=1:1:vertice_per_facets
            old = facets(line,column);
            face_(line,column) = ic(old);
        end
    end
else
    face_=[];
end


% Alternative, slow
% node_ = vertices(ia,:);
% u=unique(cells);
% for k=1:1:length(u)
%     idx = find(cells==u(k));
%     elem_(idx) = ic(u(k));
% end

% Alternative, slow
% number_vertices = length(unique_cantor);
% node_ = zeros(number_vertices,3);
% for k_vertex = 1:1:number_vertices
%     idx = find(cantor_value == unique_cantor(k_vertex));
%     node_(k_vertex,:) = vertices(idx(1),1:3);
%     for kk=1:1:length(idx)
%         elem_(cells==idx(kk))=k_vertex;
%     end
% end

    function k3 = pairing_cantor_function(k1,k2)
        k3 = 0.5.*(k1+k2).*(k1+k2+1)+k2;
    end

end

