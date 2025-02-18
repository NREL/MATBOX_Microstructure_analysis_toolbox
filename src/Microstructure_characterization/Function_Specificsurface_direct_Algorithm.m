function [Postion_surface,surface,LOC_surface] = Function_Specificsurface_direct_Algorithm(binary_phase)

% If you just want to have the surface value, you can use this one liner
% in 2D: sum(sum( abs(diff(BW,1,1)) )) + sum(sum( abs(diff(BW,1,2)) ))
% In 3D: sum(sum(sum( abs(diff(BW,1,1)) ))) + sum(sum(sum( abs(diff(BW,1,2)) ))) + sum(sum(sum( abs(diff(BW,1,3)) )))

% Domain size
Current_Domain_size = size(binary_phase);

if length(Current_Domain_size)==3
    % Position of the faces, Initialization
    % One step is equal to half voxel length
    Postion_surface.direction(1).value=zeros(2*Current_Domain_size(1)-1,2);
    Postion_surface.direction(2).value=zeros(2*Current_Domain_size(2)-1,2);
    Postion_surface.direction(3).value=zeros(2*Current_Domain_size(3)-1,2);
    % Position (axis)
    Postion_surface.direction(1).value(:,1)=(0.5:0.5:(0.5*2*Current_Domain_size(1)-0.5))';
    Postion_surface.direction(2).value(:,1)=(0.5:0.5:(0.5*2*Current_Domain_size(2)-0.5))';
    Postion_surface.direction(3).value(:,1)=(0.5:0.5:(0.5*2*Current_Domain_size(3)-0.5))';

    % The xor(A,B) function (which is the matlab built-in function for the
    % logical exclusive-OR) allows to calculate the boudnary surface easily:
    % example: the phase is [0 0 1 1 0 0 1 0]
    % then :              A=[0 0 1 1 0 0 1]
    %                     B=  [0 1 1 0 0 1 0]
    % xor(A,B)=xor( [0 0 1 1 0 0 1],
    %               [0 1 1 0 0 1 0] );
    %         =     [0 1 0 1 0 1 1]
    % suface=sum(xor(A,B))=4

    % Detection of the faces normal to the direction 1
    Face_normal_1 = xor( binary_phase(1:Current_Domain_size(1)-1,:,:) , binary_phase(2:Current_Domain_size(1),:,:));
    % Summation of all these faces
    sum1 = sum(sum(sum( Face_normal_1 )));
    % Position of all these faces
    for position_1=1:1:Current_Domain_size(1)-1
        Postion_surface.direction(1).value(2*position_1,2)= sum(sum(Face_normal_1(position_1,:,:)==1));
    end
    for position_2=1:1:Current_Domain_size(2)
        Postion_surface.direction(2).value(2*position_2-1,2)= sum(sum(Face_normal_1(:,position_2,:)==1));
    end
    for position_3=1:1:Current_Domain_size(3)
        Postion_surface.direction(3).value(2*position_3-1,2)= sum(sum(Face_normal_1(:,:,position_3)==1));
    end

    % Detection of the faces normal to the direction 2
    Face_normal_2 = xor( binary_phase(:,1:Current_Domain_size(2)-1,:) , binary_phase(:,2:Current_Domain_size(2),:));
    % Summation of all these faces
    sum2 = sum(sum(sum( Face_normal_2 )));
    % Position of all these faces
    for position_1=1:1:Current_Domain_size(1)
        Postion_surface.direction(1).value(2*position_1-1,2)= Postion_surface.direction(1).value(2*position_1-1,2) + sum(sum(Face_normal_2(position_1,:,:)==1));
    end
    for position_2=1:1:Current_Domain_size(2)-1
        Postion_surface.direction(2).value(2*position_2,2)= Postion_surface.direction(2).value(2*position_2,2) + sum(sum(Face_normal_2(:,position_2,:)==1));
    end
    for position_3=1:1:Current_Domain_size(3)
        Postion_surface.direction(3).value(2*position_3-1,2)= Postion_surface.direction(3).value(2*position_3-1,2) + sum(sum(Face_normal_2(:,:,position_3)==1));
    end

    % Detection of the faces normal to the direction 3
    Face_normal_3 = xor( binary_phase(:,:,1:Current_Domain_size(3)-1) , binary_phase(:,:,2:Current_Domain_size(3)));
    % Summation of all these faces
    sum3 = sum(sum(sum( Face_normal_3 )));
    % Position of all these faces
    for position_1=1:1:Current_Domain_size(1)
        Postion_surface.direction(1).value(2*position_1-1,2)= Postion_surface.direction(1).value(2*position_1-1,2) + sum(sum(Face_normal_3(position_1,:,:)==1));
    end
    for position_2=1:1:Current_Domain_size(2)
        Postion_surface.direction(2).value(2*position_2-1,2)= Postion_surface.direction(2).value(2*position_2-1,2) + sum(sum(Face_normal_3(:,position_2,:)==1));
    end
    for position_3=1:1:Current_Domain_size(3)-1
        Postion_surface.direction(3).value(2*position_3,2)= Postion_surface.direction(3).value(2*position_3,2) + sum(sum(Face_normal_3(:,:,position_3)==1));
    end

    % Summation of all the faces
    surface=sum1+sum2+sum3;

    % Location (for fractal dimension)
    LOC_surface=zeros(size(binary_phase));
    LOC_surface(1:Current_Domain_size(1)-1,:,:) = Face_normal_1;
    LOC_surface(:,1:Current_Domain_size(2)-1,:) = LOC_surface(:,1:Current_Domain_size(2)-1,:) + Face_normal_2;
    LOC_surface(:,:,1:Current_Domain_size(3)-1) = LOC_surface(:,:,1:Current_Domain_size(3)-1) + Face_normal_3;
    LOC_surface(LOC_surface~=0)=1;

else
    % Position of the faces, Initialization
    % One step is equal to half voxel length
    Postion_surface.direction(1).value=zeros(2*Current_Domain_size(1)-1,2);
    Postion_surface.direction(2).value=zeros(2*Current_Domain_size(2)-1,2);
    % Position (axis)
    Postion_surface.direction(1).value(:,1)=(0.5:0.5:(0.5*2*Current_Domain_size(1)-0.5))';
    Postion_surface.direction(2).value(:,1)=(0.5:0.5:(0.5*2*Current_Domain_size(2)-0.5))';

    % The xor(A,B) function (which is the matlab built-in function for the
    % logical exclusive-OR) allows to calculate the boudnary surface easily:
    % example: the phase is [0 0 1 1 0 0 1 0]
    % then :              A=[0 0 1 1 0 0 1]
    %                     B=  [0 1 1 0 0 1 0]
    % xor(A,B)=xor( [0 0 1 1 0 0 1],
    %               [0 1 1 0 0 1 0] );
    %         =     [0 1 0 1 0 1 1]
    % suface=sum(xor(A,B))=4

    % Detection of the faces normal to the direction 1
    Face_normal_1 = xor( binary_phase(1:Current_Domain_size(1)-1,:) , binary_phase(2:Current_Domain_size(1),:));
    % Summation of all these faces
    sum1 = sum(sum(sum( Face_normal_1 )));
    % Position of all these faces
    for position_1=1:1:Current_Domain_size(1)-1
        Postion_surface.direction(1).value(2*position_1,2)= sum(sum(Face_normal_1(position_1,:)==1));
    end
    for position_2=1:1:Current_Domain_size(2)
        Postion_surface.direction(2).value(2*position_2-1,2)= sum(sum(Face_normal_1(:,position_2)==1));
    end

    % Detection of the faces normal to the direction 2
    Face_normal_2 = xor( binary_phase(:,1:Current_Domain_size(2)-1) , binary_phase(:,2:Current_Domain_size(2)));
    % Summation of all these faces
    sum2 = sum(sum(sum( Face_normal_2 )));
    % Position of all these faces
    for position_1=1:1:Current_Domain_size(1)
        Postion_surface.direction(1).value(2*position_1-1,2)= Postion_surface.direction(1).value(2*position_1-1,2) + sum(sum(Face_normal_2(position_1,:)==1));
    end
    for position_2=1:1:Current_Domain_size(2)-1
        Postion_surface.direction(2).value(2*position_2,2)= Postion_surface.direction(2).value(2*position_2,2) + sum(sum(Face_normal_2(:,position_2)==1));
    end

    % Summation of all the faces
    surface=sum1+sum2;

    % Location (for fractal dimension)
    LOC_surface=zeros(size(binary_phase));
    LOC_surface(1:Current_Domain_size(1)-1,:) = Face_normal_1;
    LOC_surface(:,1:Current_Domain_size(2)-1) = LOC_surface(:,1:Current_Domain_size(2)-1) + Face_normal_2;
    LOC_surface(LOC_surface~=0)=1;

end

end