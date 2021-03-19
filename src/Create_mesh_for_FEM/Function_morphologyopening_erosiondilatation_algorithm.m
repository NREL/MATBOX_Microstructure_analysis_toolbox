function [BW] = Function_morphologyopening_erosiondilatation_algorithm(BW,erosion_distance,tolerance,dilatation_distance,distance_method)
% Morphology opening (erosion and dilation) to simplify binary image surface
% BW = Dilation( Erosion (BW) )

% Distance map
distance_map = bwdist(~BW,distance_method);
% Erosion
BW(distance_map<=erosion_distance+tolerance)=0;
% Distance map
distance_map = bwdist(BW,distance_method);
% Dilatation
BW(distance_map<=dilatation_distance+tolerance)=1;

end

