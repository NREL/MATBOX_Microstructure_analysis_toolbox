function Phase_clean_microstructure = Function_clean_microstructure_onephase(Phase_microstructure,code1,code2)


%% ALGORITHM

% Domain size
Domain_size=size(Phase_microstructure);

% Modified phase microstructure used in this algorithm
Modified_phase_microstructure = Phase_microstructure;
% Number of voxel
voxel_number=numel(Modified_phase_microstructure);

% Initialisation
Phase_clean_microstructure=zeros(Domain_size(1),Domain_size(2),Domain_size(3));

% Create a binary microstructure : 1 = current analysed phase, 0 = complementay phase
% Initialization
binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
% Binary phase of the connected microstructure
binary_phase(Modified_phase_microstructure == code1) = 1;
% Calculate
[clean_binary_phase,~] = Function_clean_microstructure_algorithm(binary_phase);
% Attribuate
Phase_clean_microstructure(clean_binary_phase==1)=code1;
Phase_clean_microstructure(clean_binary_phase~=1)=code2;


end