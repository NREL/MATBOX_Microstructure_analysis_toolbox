function [E, W, N, S, U, L] = Neighbors (Nx, Ny, Nz)
% inputs: Nx, Ny, Nz    and     outputs: E, W, N, S, U, L
% 'Nx, Ny, Nz' are voxel dimensions of the 3D domain
% 'E, W, N, S, U, L' identify indices of E(ast), W(est), N(orth), S(outh), 
% U(pper) and L(ower) neighbors respectively
%__________________________________________________________________________
% Programmed by     : Aashutosh Mistry
% Created on        : Jul 05, 2020
% Modifications     : (none)
%__________________________________________________________________________

E = zeros(Nx,Ny,Nz);                        % east neighbor
W = zeros(Nx,Ny,Nz);                        % west neighbor
N = zeros(Nx,Ny,Nz);                        % north neighbor
S = zeros(Nx,Ny,Nz);                        % south neighbor
U = zeros(Nx,Ny,Nz);                        % upper neighbor
L = zeros(Nx,Ny,Nz);                        % lower neighbor

for k=1:Nz
    for j=1:Ny
        for i=1:Nx            
            %east neighbor
            if i<Nx
                E(i,j,k) = Index1D((i+1),j,k, Nx,Ny,Nz);
            else
                E(i,j,k) = Index1D(1,j,k, Nx,Ny,Nz);
            end
            
            %west neighbor
            if i>1
                W(i,j,k) = Index1D((i-1),j,k, Nx,Ny,Nz);
            else
                W(i,j,k) = Index1D(Nx,j,k, Nx,Ny,Nz);
            end
            
            %north neighbor
            if j<Ny
                N(i,j,k) = Index1D(i,(j+1),k, Nx,Ny,Nz);
            else
                N(i,j,k) = Index1D(i,1,k, Nx,Ny,Nz);
            end
            
            %south neighbor
            if j>1
                S(i,j,k) = Index1D(i,(j-1),k, Nx,Ny,Nz);
            else
                S(i,j,k) = Index1D(i,Ny,k, Nx,Ny,Nz);
            end
            
            %upper neighbor
            if k<Nz
                U(i,j,k) = Index1D(i,j,(k+1), Nx,Ny,Nz);
            else
                U(i,j,k) = Index1D(i,j,1, Nx,Ny,Nz);
            end
            
            %lower neighbor
            if k>1
                L(i,j,k) = Index1D(i,j,(k-1), Nx,Ny,Nz);
            else
                L(i,j,k) = Index1D(i,j,Nz, Nx,Ny,Nz);
            end            
        end
    end
end


end