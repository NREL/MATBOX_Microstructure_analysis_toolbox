function [Phase_microstructure] = Function_erosion_dilatation(Phase_microstructure,code,complementary_code)
   
    % EROSTION AND DILATATION
    % Opening(Image) = Dilation( Erosion (Image) )
    
    % Domain size
    Domain_size=size(Phase_microstructure);
    
    % Create binary matrix
    binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
    binary_phase(Phase_microstructure == code) = 1;    
        
    % Distance map
    distance_map = bwdist(~binary_phase,'chessboard');
    % Erosion
    binary_phase(distance_map<=1)=0;
    % Distance map
    distance_map = bwdist(binary_phase,'chessboard');    
     % Dilatation
    binary_phase(distance_map<=1)=1;   
    
    % Reassign
    Phase_microstructure=zeros(Domain_size(1),Domain_size(2),Domain_size(3))+complementary_code;
    Phase_microstructure(binary_phase==1)=code;    
    
end

