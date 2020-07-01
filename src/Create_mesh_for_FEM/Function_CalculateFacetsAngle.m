function [Angle_faces, Angle_verteces, Area_faces, Normal_faces, center_facets] = Function_CalculateFacetsAngle(all_verteces,all_faces, all_cell, cell_centroid)

% Number of verteces and faces
n_verteces=length(all_verteces);
n_faces=length(all_faces);

% Initialisation: angles
Angle_faces=zeros(n_faces,3);
count_verteces=zeros(n_verteces,1); % Count the times a vertece is called
Angle_verteces=zeros(n_verteces,20,2);
Area_faces=zeros(n_faces,1); % Initialisation: area
Normal_faces=zeros(n_faces,3); % Initialisation: normal
center_facets = zeros(n_faces,3); % Initialisation: center

% Loop over all faces
for k=1:1:round(n_faces/1)
    index_point_A=all_faces(k,1);
    index_point_B=all_faces(k,2);
    index_point_C=all_faces(k,3);
      
    % Coordinates of the three points
    A = all_verteces(index_point_A,:);
    B = all_verteces(index_point_B,:);
    C = all_verteces(index_point_C,:);
      
    % Coordinates of the vectors
    AB=B-A;
    AC=C-A;
    x1=AB(1); x2=AB(2); x3=AB(3);
    y1=AC(1); y2=AC(2); y3=AC(3);
    
    % Normal of the face
    nx=(x2*y3)-(x3*y2);
    ny=(x3*y1)-(x1*y3);
    nz=(x1*y2)-(x2*y1);

    % Calculation of the triange area
    % https://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle
    Area_faces(k,1)=1/2 * (nx^2 + ny^2 + nz^2)^(0.5);
    
    % Calculation of the angle (origin A)
    % https://www.mathworks.com/matlabcentral/answers/90668-how-to-calculate-a-angle-between-two-vectors-in-3d
    %angle_originA = rad2deg(atan2(norm(cross(AB,AC)),dot(AB,AC)));
    % https://www.cs.berkeley.edu/~wkahan/Mindless.pdf), section 12 "Mangled Angles."
    angle_originA= rad2deg( 2*atan2(  norm( AB*norm(AC) - norm(AB)*AC) , norm( AB*norm(AC) + norm(AB)*AC) ) );

    % Calculation of the angle (origin B)
    BA=A-B;
    BC=C-B;
    %angle_originB = rad2deg(atan2(norm(cross(BA,BC)),dot(BA,BC)));
    angle_originB = rad2deg( 2*atan2(  norm( BA*norm(BC) - norm(BA)*BC) , norm( BA*norm(BC) + norm(BA)*BC) ) );
    
    % Calculation of the angle (origin C)
    angle_originC = 180-angle_originA-angle_originB;
    %CA=A-C;
    %CB=B-C;    
    %angle_originC = rad2deg( 2*atan2(  norm( CA*norm(CB) - norm(CA)*CB) , norm( CA*norm(CB) + norm(CA)*CB) ) );
       
    % All angle, per face
    Angle_faces(k,1)=angle_originA;
    Angle_faces(k,2)=angle_originB;
    Angle_faces(k,3)=angle_originC;

end

% plane=surfplane(all_verteces,all_faces);
% norm_ = abs(plane(:,1))+abs(plane(:,2))+abs(plane(:,3));
% Normal_faces(:,1)=plane(:,1)./norm_;
% Normal_faces(:,2)=plane(:,2)./norm_;
% Normal_faces(:,3)=plane(:,3)./norm_;

%  snorm=surfacenorm(node,face)
%     or
%  snorm=surfacenorm(node,face,'Normalize',0)
%  compute the normal vectors for a triangular surface
%  input:
%    node: a list of node coordinates (nn x 3)
%    face: a surface mesh triangle list (ne x 3)
%    opt: a list of optional parameters, currently surfacenorm supports:
%         'Normalize': [1|0] if set to 1, the normal vectors will be 
%                            unitary (default)
%  output:
%    snorm: output surface normal vector at each face

%[newnode,newface]=surfreorient(all_verteces,all_cell);
%Normal_faces=surfacenorm(newnode,newface);


end

