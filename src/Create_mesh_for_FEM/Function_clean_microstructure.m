function Phase_clean_microstructure = Function_clean_microstructure(Phase_microstructure,Domain_size)


%% GLOBAl VARIABLES
global INFO;
global OPTIONS;


%%
%% START DATE
%%
%%



if OPTIONS.date==1;
    date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
    if OPTIONS.display_textresults==1;
        fprintf ('   Started the %s\n\n',date_start);
    end
end



%% ALGORITHM

% CPU and stopwatch time - start
time_cpu_start = cputime;
tic;

% Modified phase microstructure used in this algorithm
Modified_phase_microstructure = Phase_microstructure;

% Number of phase
number_phase = length(INFO.phase);

% Number of voxel
voxel_number=numel(Modified_phase_microstructure);

% Initialisation
Phase_clean_microstructure=zeros(Domain_size(1),Domain_size(2),Domain_size(3));

% Calculation
for current_phase=1:1:number_phase
    
    if current_phase>1
       % We will reallocate the voxels removed from the previous phase to the next phases
       % These 'removed voxels' corresponds only to the voxel that have
       % been removed from the line and node 'only' connection detection.
       % The voxels that have been removed due to a potential loss of
       % connectivity are not concerned here.
       
       % Values of the previous phase
       number_previous_phase = current_phase-1;
       code_previous_phase=zeros(1,number_previous_phase);
       for previous_phase=current_phase-1:-1:1
           code_previous_phase(1,previous_phase) = INFO.phase(previous_phase).code;
       end
              
       number_removed_voxel = length(index_removed_voxels);
       if number_removed_voxel>0
           % Get back coordinates of all these voxels
           [IX,IY,IZ]= ind2sub(Domain_size,index_removed_voxels);
           for current_voxel=1:1:number_removed_voxel
               % Coordinate of the current voxel
               x_=IX(current_voxel); y_=IY(current_voxel); z_=IZ(current_voxel);
               % Get surroundin values (only face-face connection)
               surrounding_value = [];
               if x_>1
                   surrounding_value =[surrounding_value Modified_phase_microstructure(x_-1,y_,z_)];
               end
               if y_>1
                   surrounding_value =[surrounding_value Modified_phase_microstructure(x_,y_-1,z_)];
               end               
                if z_>1
                   surrounding_value =[surrounding_value Modified_phase_microstructure(x_,y_,z_-1)];
               end              
               if x_<Domain_size(1)
                   surrounding_value =[surrounding_value Modified_phase_microstructure(x_+1,y_,z_)];
               end               
               if y_<Domain_size(2)
                   surrounding_value =[surrounding_value Modified_phase_microstructure(x_,y_+1,z_)];
               end      
               if z_<Domain_size(3)
                   surrounding_value =[surrounding_value Modified_phase_microstructure(x_,y_,z_+1)];
               end                     
               surrounding_value=real(surrounding_value); % Only phase information
               unique_surrounding = unique(surrounding_value); % Different surrounding phase
               n_unique = length(unique_surrounding); % Number of dufferent phase around
               sum_surrounding = zeros(n_unique,2); % Initialise
               for k=1:1:n_unique % Counting phase voxel around
                   sum_surrounding(k,1)=unique_surrounding(k);
                   sum_surrounding(k,2)=sum(surrounding_value==unique_surrounding(k));
               end
               sum_surrounding = sortrows(sum_surrounding,-2); % sort by decreasing order
               % Remove previous phase
               for previous_phase=1:1:number_previous_phase
                   index_previous = find(sum_surrounding(:,1)==code_previous_phase(1,previous_phase));
                   sum_surrounding(index_previous,:)=[];
               end
               % Re-attribute the voxel
               if  ~isempty(sum_surrounding)
                   Modified_phase_microstructure(x_,y_,z_)=sum_surrounding(1,1)+1i;
               end
           end
       end
    end
    
    % Code of the phase
    code_tmp = INFO.phase(current_phase).code;
    % Create a binary microstructure : 1 = current analysed phase, 0 = complementay phase
    % Initialization
    binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
    % Binary phase of the connected microstructure
    binary_phase(Modified_phase_microstructure == (code_tmp+1i)) = 1;
    % Calculate
    [clean_binary_phase,index_removed_voxels] = Function_clean_microstructure_algorithm(binary_phase);
    % Attribuate
    Phase_clean_microstructure(clean_binary_phase==1)=code_tmp+1i;
end

% CPU and stopwatch time - end
time_cpu_elapsed = cputime-time_cpu_start;
time_stopwatch_elapsed = toc;
Time_measure = [voxel_number time_cpu_elapsed time_stopwatch_elapsed];




end