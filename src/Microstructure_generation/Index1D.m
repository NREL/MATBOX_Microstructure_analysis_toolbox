function [P] = Index1D (i,j,k, Nx,Ny,Nz)
% inputs: i,j,k, Nx,Ny,Nz    and     outputs: P
% 'i, j, k' are 3D voxel coordinates
% 'Nx, Ny, Nz' are voxel dimensions of the 3D domain
% 'P' is corresponding 1D coordinate
%__________________________________________________________________________
% Programmed by     : Aashutosh Mistry
% Created on        : Jul 05, 2020
% Modifications     : (none)
%__________________________________________________________________________

P = i + ((j-1)*Nx) + ((k-1)*Nx*Ny);

end