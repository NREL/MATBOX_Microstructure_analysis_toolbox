function [Connectivity_matrix, Interface_label] = Function_connectivitymatrix(M, background, interface_complemetaryvolume, interface_anotherlabel)

Domain_size = size(M);
dimension = length(Domain_size);
Interface_label = zeros(Domain_size);

Ids = unique(M); % Get all ids
Ids(Ids==background) = []; % Remove background
n_ids = length(Ids);

% Connectivity initialisation
Connectivity_matrix = zeros(n_ids+1,n_ids+1);
% Set table axis
for k=1:1:n_ids
    id = Ids(k); % Get current id
    % Write axis
    Connectivity_matrix(k+1,1)=id;
    Connectivity_matrix(1,k+1)=id;
    % Write identity
    Connectivity_matrix(k+1,k+1)=1;
end

% Connectivity_matrix(i,j) is the number of voxels identified at the interface (if any) between the two labels
% with a face-to-face identification for the interface, then the number is equal to the number of voxel face at the interface
number_voxel_atinterface=0; % Initialize
% loop over all labels
for k=1:1:n_ids
    % Get id_
    id = Ids(k);
    % Find all voxels of the particle
    index_label=find(M==id);
    % Get all their coordinate
    [Pa1,Pa2,Pa3] = ind2sub(Domain_size,index_label);
    % Number of voxel
    number_voxel_label = length(Pa1);
    % Loop over all voxels of the phase
    for current_voxel=1:1:number_voxel_label
        % Get back voxel coordinate
        x_=Pa1(current_voxel); y_=Pa2(current_voxel); z_=Pa3(current_voxel);
        
        % Check x-
        if x_>1
            % Tested voxel
            x_tested = x_-1; y_tested = y_; z_tested = z_;
            % Value founded there
            Id_tested = M(x_tested,y_tested,z_tested);
            if Id_tested~=id
                % Interface founded
                if Id_tested==background
                    % Interface with the complementary phase
                    Interface_label(x_,y_,z_)=interface_complemetaryvolume;
                else
                    % Interface with the another particle of the phase
                    Interface_label(x_,y_,z_)=interface_anotherlabel;
                    % Update Connectivity_matrix as consequence
                    row_ = find( Connectivity_matrix(:,1)==id);
                    column_ = find( Connectivity_matrix(1,:)==Id_tested);
                    Connectivity_matrix(row_,column_)=Connectivity_matrix(row_,column_)+1;
                end
            end
        end
        
        % Check y-
        if y_>1
            % Tested voxel
            x_tested = x_; y_tested = y_-1; z_tested = z_;
            % Value founded there
            Id_tested = M(x_tested,y_tested,z_tested);
            if Id_tested~=id
                % Interface founded
                if Id_tested==background
                    % Interface with the complementary phase
                    Interface_label(x_,y_,z_)=interface_complemetaryvolume;
                else
                    % Interface with the another particle of the phase
                    Interface_label(x_,y_,z_)=interface_anotherlabel;
                    % Update Connectivity_matrix as consequence
                    row_ = find( Connectivity_matrix(:,1)==id);
                    column_ = find( Connectivity_matrix(1,:)==Id_tested);
                    Connectivity_matrix(row_,column_)=Connectivity_matrix(row_,column_)+1;
                end
            end
        end
        
        % Check z-
        if dimension==3
            if z_>1
                % Tested voxel
                x_tested = x_; y_tested = y_; z_tested = z_-1;
                % Value founded there
                Id_tested = M(x_tested,y_tested,z_tested);
                if Id_tested~=id
                    % Interface founded
                    if Id_tested==background
                        % Interface with the complementary phase
                        Interface_label(x_,y_,z_)=interface_complemetaryvolume;
                    else
                        % Interface with the another particle of the phase
                        Interface_label(x_,y_,z_)=interface_anotherlabel;
                        % Update Connectivity_matrix as consequence
                        row_ = find( Connectivity_matrix(:,1)==id);
                        column_ = find( Connectivity_matrix(1,:)==Id_tested);
                        Connectivity_matrix(row_,column_)=Connectivity_matrix(row_,column_)+1;
                    end
                end
            end
        end
        
        % Check x+
        if x_<Domain_size(1)
            % Tested voxel
            x_tested = x_+1; y_tested = y_; z_tested = z_;
            % Value founded there
            Id_tested = M(x_tested,y_tested,z_tested);
            if Id_tested~=id
                % Interface founded
                if Id_tested==background
                    % Interface with the complementary phase
                    Interface_label(x_,y_,z_)=interface_complemetaryvolume;
                else
                    % Interface with the another particle of the phase
                    Interface_label(x_,y_,z_)=interface_anotherlabel;
                    % Update Connectivity_matrix as consequence
                    row_ = find( Connectivity_matrix(:,1)==id);
                    column_ = find( Connectivity_matrix(1,:)==Id_tested);
                    Connectivity_matrix(row_,column_)=Connectivity_matrix(row_,column_)+1;
                end
            end
        end
        
        % Check y+
        if y_<Domain_size(2)
            % Tested voxel
            x_tested = x_; y_tested = y_+1; z_tested = z_;
            % Value founded there
            Id_tested = M(x_tested,y_tested,z_tested);
            if Id_tested~=id
                % Interface founded
                if Id_tested==background
                    % Interface with the complementary phase
                    Interface_label(x_,y_,z_)=interface_complemetaryvolume;
                else
                    % Interface with the another particle of the phase
                    Interface_label(x_,y_,z_)=interface_anotherlabel;
                    % Update Connectivity_matrix as consequence
                    row_ = find( Connectivity_matrix(:,1)==id);
                    column_ = find( Connectivity_matrix(1,:)==Id_tested);
                    Connectivity_matrix(row_,column_)=Connectivity_matrix(row_,column_)+1;
                end
            end
        end
        
        % Check z+
        if dimension==3
            if z_<Domain_size(3)
                % Tested voxel
                x_tested = x_; y_tested = y_; z_tested = z_+1;
                % Value founded there
                Id_tested = M(x_tested,y_tested,z_tested);
                if Id_tested~=id
                    % Interface founded
                    if Id_tested==background
                        % Interface with the complementary phase
                        Interface_label(x_,y_,z_)=interface_complemetaryvolume;
                    else
                        % Interface with the another particle of the phase
                        Interface_label(x_,y_,z_)=interface_anotherlabel;
                        % Update Connectivity_matrix as consequence
                        row_ = find( Connectivity_matrix(:,1)==id);
                        column_ = find( Connectivity_matrix(1,:)==Id_tested);
                        Connectivity_matrix(row_,column_)=Connectivity_matrix(row_,column_)+1;
                    end
                end
            end
        end
        
    end
end
% Check the matrix is symmetric
if ~issymmetric(Connectivity_matrix)
    warning('Connectivity matrix is non-symmetrical (Error)')
end

end

