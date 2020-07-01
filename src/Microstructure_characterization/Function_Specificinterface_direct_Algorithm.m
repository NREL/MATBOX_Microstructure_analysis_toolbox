function [Postion_surface,surface] = Function_Specificinterface_direct_Algorithm(binary_phase)

% Domain size
Current_Domain_size = size(binary_phase);

% Position of the faces, Initialization
% One step is equal to half voxel length
Postion_surface.direction(1).value=zeros(2*Current_Domain_size(1)-1,2);
Postion_surface.direction(2).value=zeros(2*Current_Domain_size(2)-1,2);
Postion_surface.direction(3).value=zeros(2*Current_Domain_size(3)-1,2);
% Position (axis)
Postion_surface.direction(1).value(:,1)=(0.5:0.5:(0.5*2*Current_Domain_size(1)-0.5))';
Postion_surface.direction(2).value(:,1)=(0.5:0.5:(0.5*2*Current_Domain_size(2)-0.5))';
Postion_surface.direction(3).value(:,1)=(0.5:0.5:(0.5*2*Current_Domain_size(3)-0.5))';

method_fast = true;

if method_fast

    % Detection of the faces normal to the direction 1
    Face_normal_1 = abs (binary_phase(1:Current_Domain_size(1)-1,:,:) - binary_phase(2:Current_Domain_size(1),:,:));
    Face_normal_1(Face_normal_1~=9)=0;
    Face_normal_1(Face_normal_1==9)=1;
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
    Face_normal_2 = abs( binary_phase(:,1:Current_Domain_size(2)-1,:) - binary_phase(:,2:Current_Domain_size(2),:));
    Face_normal_2(Face_normal_2~=9)=0;
    Face_normal_2(Face_normal_2==9)=1;    
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
    Face_normal_3 = abs( binary_phase(:,:,1:Current_Domain_size(3)-1) - binary_phase(:,:,2:Current_Domain_size(3)));
    Face_normal_3(Face_normal_3~=9)=0;
    Face_normal_3(Face_normal_3==9)=1;      
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
    
else
    
    % Initialisation
    sum1=0; sum2=0; sum3=0;
    
    % Detection of the faces normal to the direction 1
    for position_1=1:1:Current_Domain_size(1)-1
        for position_2=1:1:Current_Domain_size(2)
            for position_3=1:1:Current_Domain_size(3)
                voxel_A = binary_phase(position_1,position_2,position_3);
                voxel_B = binary_phase(position_1+1,position_2,position_3);
                if abs(voxel_A-voxel_B)==9
                    sum1=sum1+1;
                    Postion_surface.direction(1).value(2*position_1,2)= Postion_surface.direction(1).value(2*position_1,2)+1;
                    Postion_surface.direction(2).value(2*position_2-1,2)= Postion_surface.direction(2).value(2*position_2-1,2)+1;
                    Postion_surface.direction(3).value(2*position_3-1,2)= Postion_surface.direction(3).value(2*position_3-1,2)+1;
                end
            end
        end
    end

    % Detection of the faces normal to the direction 2
    for position_2=1:1:Current_Domain_size(2)-1
        for position_1=1:1:Current_Domain_size(1)
            for position_3=1:1:Current_Domain_size(3)
                voxel_A = binary_phase(position_1,position_2,position_3);
                voxel_B = binary_phase(position_1,position_2+1,position_3);
                if abs(voxel_A-voxel_B)==9
                    sum2=sum2+1;
                    Postion_surface.direction(1).value(2*position_1-1,2)= Postion_surface.direction(1).value(2*position_1-1,2)+1;
                    Postion_surface.direction(2).value(2*position_2,2)= Postion_surface.direction(2).value(2*position_2,2)+1;
                    Postion_surface.direction(3).value(2*position_3-1,2)= Postion_surface.direction(3).value(2*position_3-1,2)+1;
                end
            end
        end
    end

    % Detection of the faces normal to the direction 3
    for position_3=1:1:Current_Domain_size(3)-1
        for position_1=1:1:Current_Domain_size(1)
            for position_2=1:1:Current_Domain_size(2)
                voxel_A = binary_phase(position_1,position_2,position_3);
                voxel_B = binary_phase(position_1,position_2,position_3+1);
                if abs(voxel_A-voxel_B)==9
                    sum3=sum3+1;
                    Postion_surface.direction(1).value(2*position_1-1,2)= Postion_surface.direction(1).value(2*position_1-1,2)+1;
                    Postion_surface.direction(2).value(2*position_2-1,2)= Postion_surface.direction(2).value(2*position_2-1,2)+1;
                    Postion_surface.direction(3).value(2*position_3,2)= Postion_surface.direction(3).value(2*position_3,2)+1;
                end
            end
        end
    end

    % Summation of all the faces
    surface=sum1+sum2+sum3;
  
end
    
end