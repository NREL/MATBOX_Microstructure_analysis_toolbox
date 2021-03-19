function [I] = AddCBDwMorphology_step (I, w, Pid,Eid,Wid,Nid,Sid,Uid,Lid, Nsimultaneous)
% inputs: I, w, Pid,Eid,Wid,Nid,Sid,Uid,Lid     and     output: I
% 'I' is Nx.Ny.Nz array with 0 as pore, 1 as active material, 2 as CBD
% 'w' is corresponding morphology (varies between 0 and 1)
% '*id' are 1D indices for each voxel (P) and its neighbors (E,W,N,S,U,L)
% 'Nsimultaneous' many voxels are added in one execution (exact number may
% be smaller)
%__________________________________________________________________________
% Programmed by     : Aashutosh Mistry
% Created on        : Jul 05, 2020
% Modifications     : Aug 29, 2020
%                       -   candidate voxels are void voxels having
%                           at least one solid neighbor
%__________________________________________________________________________

[Nx,Ny,Nz] = size(I);                       % voxels in different directions
N = Nx*Ny*Nz;                               % total voxel counts

kmin = 2; kmax = Nz-1; normalize_energy = 6;
if Nz==1 % 2D case
    kmin = 1; kmax = 1; normalize_energy = 4;
end

% flattening 'I' matrix into 1D array _____________________________________
I1D = zeros(N,1);

% for k=1:Nz
%     for j=1:Ny
%         for i=1:Nx
%             I1D(Pid(i,j,k)) = I(i,j,k);
%         end
%     end
% end
I1D(Pid) = I;

% candidate voxels : void voxels with at least one solid neighbor _________
Ncandidate=0;

for k=kmin:kmax
    for j=2:Ny-1
        for i=2:Nx-1
            if I(i,j,k)==0
                % reinitialize
                ngbr1 = 0;                      % active material neighbors
                ngbr2 = 0;                      % CBD neighbors

                % checking East neighbor
                if I1D(Eid(i,j,k))==1
                    ngbr1 = ngbr1 + 1;
                elseif I1D(Eid(i,j,k))==2
                    ngbr2 = ngbr2 + 1;
                end

                % checking West neighbor
                if I1D(Wid(i,j,k))==1
                    ngbr1 = ngbr1 + 1;
                elseif I1D(Wid(i,j,k))==2
                    ngbr2 = ngbr2 + 1;
                end

                % checking North neighbor
                if I1D(Nid(i,j,k))==1
                    ngbr1 = ngbr1 + 1;
                elseif I1D(Nid(i,j,k))==2
                    ngbr2 = ngbr2 + 1;
                end

                % checking South neighbor
                if I1D(Sid(i,j,k))==1
                    ngbr1 = ngbr1 + 1;
                elseif I1D(Sid(i,j,k))==2
                    ngbr2 = ngbr2 + 1;
                end

                % checking Upper neighbor
                if I1D(Uid(i,j,k))==1
                    ngbr1 = ngbr1 + 1;
                elseif I1D(Uid(i,j,k))==2
                    ngbr2 = ngbr2 + 1;
                end

                % checking Lower neighbor
                if I1D(Lid(i,j,k))==1
                    ngbr1 = ngbr1 + 1;
                elseif I1D(Lid(i,j,k))==2
                    ngbr2 = ngbr2 + 1;
                end

                % if there are solid neighbors, qualifies as a candidate voxel
                if (ngbr1+ngbr2)>0
                    Ncandidate = Ncandidate + 1;
                    Qid(Ncandidate) = Pid(i,j,k);                                       % index of candidate voxel
                    Energy(Ncandidate) = (ngbr2*w/normalize_energy) + (ngbr1*(1-w)/normalize_energy);                 % eneregy gain with deposition
                    if Ncandidate>1
                        cumEnergy(Ncandidate) = cumEnergy(Ncandidate-1) + Energy(Ncandidate);
                    else
                        cumEnergy(Ncandidate) = Energy(Ncandidate);
                    end
                end
            end
        end
    end
end

% selecting deposition voxels based on energy map _________________________
% Selected_energy = zeros(Nsimultaneous,1) - 1;
for m=1:Nsimultaneous
    target_cumEnergy = cumEnergy(Ncandidate)*rand();
    
    for n=2:Ncandidate
        if ((target_cumEnergy>cumEnergy(n-1)) && (target_cumEnergy<=cumEnergy(n)))
            I1D(Qid(n)) = 2;
            % Selected_energy(m,1) = Energy(n);
            break;
        end
    end
end
% Selected_energy( Selected_energy==-1) = []; % Remove non assigned value
% mean(Selected_energy) / mean(Energy) % > 1 means selected voxel on average have higher than average energy deposition


% converting back to 3D 'I' matrix ________________________________________
% for k=1:Nz
%     for j=1:Ny
%         for i=1:Nx
%             I(i,j,k) = I1D(Pid(i,j,k));
%         end
%     end
% end
I = I1D(Pid);

end