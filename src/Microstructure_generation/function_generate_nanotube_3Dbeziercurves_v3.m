function [tube_phase, tube_id, tube_skeleton, tubes, outcome] = function_generate_nanotube_3Dbeziercurves_v3(domain_size,Bezier_pars,sav,stopcond,cropp,p)

%%
% x,y,z orientation, no variation, on edges: isotropy, random tube 
% x,y,z orientation, no variation, center: isotropy, random tube 

% curved, no variation, on edges:
% 0 elevation, x,y isotropy, random tube


%%

visualize_seeds = false; % True for debugging
tube_can_overlap_with_mask = false; % By default

Qprecision = 10; % Set high value to avoid tube non-contiguity

check_isotropy = false;
if check_isotropy
    k_iso = 0;
    tau_tube = [];
    aniso_tube = [];
    tau_pore = [];
    aniso_pore = [];
end

%% THIRD-PARTY LICENSE
% I am using a slighly modified function "hobbysplines" to create smooth 3D bezier curves from:

% Will Robertson (2023). Smooth 3D bezier curves with implicit control points
% https://www.mathworks.com/matlabcentral/fileexchange/42302-smooth-3d-bezier-curves-with-implicit-control-points
% MATLAB Central File Exchange. Retrieved March 29, 2023.
% License is available within custom_hobbysplines.m

% Everything else is NREL code covered with the MATBOX licence

%% INITIALIZATION
tube_phase = zeros(domain_size);
tube_id = zeros(domain_size);
tube_skeleton = zeros(domain_size);
if isempty(p.maskused)
    p.maskused = ones(domain_size);
end

%% SET TUBE SEEDS
% =1 possible starting location
% =0 forbidden starting location

% Initialization
seeds = tube_phase+1;
used_seeds = zeros(domain_size);

if p.seeds_onlyatedges % Locate seeds at domain's edge
    seeds(2:end-1,2:end-1,2:end-1)=0;
    if sum(cropp.distfromedge)>0
        l1 = round(domain_size(1)*cropp.distfromedge(1));
        l2 = round(domain_size(2)*cropp.distfromedge(2));
        l3 = round(domain_size(3)*cropp.distfromedge(3));

        seeds(:,1:l2,1:l3)=0; seeds(:,1:l2,end-l3+1:end)=0;
        seeds(:,end-l2+1:end,1:l3)=0; seeds(:,end-l2+1:end,end-l3+1:end)=0;

        seeds(1:l1,:,1:l3)=0; seeds(1:l1,:,end-l3+1:end)=0;
        seeds(end-l1+1:end,:,1:l3)=0; seeds(end-l1+1:end,:,end-l3+1:end)=0;

        seeds(1:l1,1:l2,:)=0; seeds(1:l1,end-l2+1:end,:)=0;
        seeds(end-l1+1:end,1:l2,:)=0; seeds(end-l1+1:end,end-l2+1:end,:)=0;
    end
end

seeds = seeds .* p.maskused; % Starting location can only be within mask

idseeds = find(seeds); % Index of all initial starting locations

if visualize_seeds
    Microstructure_basic_visualization_interface(seeds)
end

%% ALGORITHM

if p.seeds_onlyatedges && p.force_seedequidistribution
    seed_xyz = [0 0 0];
end

nvoxel = numel(tube_phase);
% Initalize wallclock time and stoping conditions
tStart_total = tic;
tStart_tube = tic;
time_currenttube = 0;
k_tube = 0;
volume_fraction = 0;
outcome = 'Success'; % By default
while k_tube < stopcond.n_tube && volume_fraction < stopcond.target_volumefraction
    t_total = toc(tStart_total);

    % Check time related stoping condition
    if t_total > stopcond.totaltime
        outcome = 'Exceed total time';
        return;
    end
    if time_currenttube > stopcond.pertube
        outcome = 'Exceed tube time';
        return;
    end

    % Select radius, tube length, initial orientation within a distribution
    radius = pickfromdistriubtion(p.maxdeviation_radius,p.equiprobability_radius,p.mean_radius,p.cdf_radius,p.x_radius);
    if p.finitelength
        tubelength = pickfromdistriubtion(p.maxdeviation_length,p.equiprobability_length,p.mean_length,p.cdf_length,p.x_length);
    end
    if  strcmp(p.initial_orientation_choice,'Angle')
        azimuth = pickfromdistriubtion(p.maxdeviation_initialazimuth,p.equiprobability_initialazimuth,p.mean_initialazimuth,p.cdf_initialazimuth,p.x_initialazimuth);
        elevation = pickfromdistriubtion(p.maxdeviation_initialelevation,p.equiprobability_initialelevation,p.mean_initialelevation,p.cdf_initialelevation,p.x_initialelevation);
    else
        Ux = pickfromdistriubtion(p.maxdeviation_initialUx,p.equiprobability_initialUx,p.mean_initialUx,p.cdf_initialUx,p.x_initialUx);
        Uy = pickfromdistriubtion(p.maxdeviation_initialUy,p.equiprobability_initialUy,p.mean_initialUy,p.cdf_initialUy,p.x_initialUy);
        Uz = pickfromdistriubtion(p.maxdeviation_initialUz,p.equiprobability_initialUz,p.mean_initialUz,p.cdf_initialUz,p.x_initialUz);
        %r = rand(3,1);
        %u = r/norm(r);
        %[azimuth,elevation,~] = cart2sph(u(1),u(2),u(3));
        [azimuth,elevation,~] = cart2sph(Ux,Uy,Uz);
        azimuth =rad2deg(azimuth);
        elevation =rad2deg(elevation);        
    end

    % Select seed
    if p.seeds_onlyatedges && k_tube>0 && p.force_seedequidistribution
        [IX,IY,IZ]=ind2sub(domain_size,idseeds);
        idx = find(double(IX==1) + double(IX==domain_size(1)));
        idy = find(double(IY==1) + double(IY==domain_size(2)));
        idz = find(double(IZ==1) + double(IZ==domain_size(3)));
        
        ind = find(seed_xyz==min(seed_xyz));
        ind=ind(1); % In case of equality
        if ind==1
            idseed = max([1 round(numel(idx)*rand)]);
            x_previous=IX(idx(idseed)); y_previous=IY(idx(idseed)); z_previous=IZ(idx(idseed));
        elseif ind==2
            idseed = max([1 round(numel(idy)*rand)]);
            x_previous=IX(idy(idseed)); y_previous=IY(idy(idseed)); z_previous=IZ(idy(idseed));
        elseif ind==3
            idseed = max([1 round(numel(idz)*rand)]);
            x_previous=IX(idz(idseed)); y_previous=IY(idz(idseed)); z_previous=IZ(idz(idseed));
        end

    else
        idseed = max([1 round(numel(idseeds)*rand)]);
        [x_previous, y_previous, z_previous] = ind2sub(domain_size,idseeds(idseed));
    end


    if p.seeds_onlyatedges && p.force_seedequidistribution
        if x_previous==1 || x_previous==domain_size(1)
            seed_candidate = [1 0 0];
        elseif y_previous==1 || y_previous==domain_size(2)
            seed_candidate = [0 1 0];
        elseif z_previous==1 || z_previous==domain_size(3)
            seed_candidate = [0 0 1];
        else
            warning('incorrect seed location');
        end
    end

    % Seed coordinate
    x_seed = x_previous;
    y_seed = y_previous;
    z_seed = z_previous;

    straight_tube = false;
    if strcmp(p.change_orientation_choice,'Angle')
        if p.maxdeviation_azimuthchange==0 && p.mean_azimuthchange ==0 && p.maxdeviation_elevationchange==0 && p.mean_elevationchange==0
            straight_tube = true;
        end
    else
        if p.maxdeviation_Uxchange==0 && p.mean_Uxchange ==0 && p.maxdeviation_Uychange==0 && p.mean_Uychange ==0 && p.maxdeviation_Uzchange==0 && p.mean_Uzchange ==0
            straight_tube = true;
        end
    end

    if straight_tube % Straight tube, simplified code
        % We need only two points: start and end
        if p.finitelength
            distancebetweenangleupdate = tubelength;
        else
            distancebetweenangleupdate = p.min_length;
        end
        sense = +1;
        [x_step,y_step,z_step] = sph2cart2(deg2rad(azimuth),deg2rad(elevation),distancebetweenangleupdate);
        x_new =  round(x_previous+x_step);
        y_new =  round(y_previous+y_step);
        z_new =  round(z_previous+z_step);
        if x_new<1 || y_new<1 || z_new<1 || x_new>domain_size(1) || y_new>domain_size(2) || z_new>domain_size(3) % Out of bounds
            x_new =  round(x_previous-x_step);
            y_new =  round(y_previous-y_step);
            z_new =  round(z_previous-z_step);
            sense = -1;
            if x_new<1 || y_new<1 || z_new<1 || x_new>domain_size(1) || y_new>domain_size(2) || z_new>domain_size(3) % Out of bounds
                time_currenttube = toc(tStart_tube);
                continue % Pick another seed and another orientation
            end
        end

        % If tube are infinite in length and end point is not at domain's edge, push it there
        if ~p.finitelength && (x_new~=1 && y_new~=1 && z_new~=1 && x_new~=domain_size(1) && y_new~=domain_size(2) && z_new~=domain_size(3)) % Tube didnt reach edge
            slopes = [(x_new-x_previous) (y_new-y_previous) (z_new-z_previous)]/distancebetweenangleupdate;
            backward = [(1-x_previous) (1-y_previous) (1-z_previous)]./slopes;
            forward = [(domain_size(1)-x_previous) (domain_size(2)-y_previous) (domain_size(3)-z_previous)]./slopes;
            backward(backward<0)=0;
            forward(forward<0)=0;
            t = backward+forward;
            idx = t==min(t);
            distancebetweenangleupdate = t(idx)+1;
            distancebetweenangleupdate = distancebetweenangleupdate(1); % In case of equality
            [x_step,y_step,z_step] = sph2cart2(deg2rad(azimuth),deg2rad(elevation),distancebetweenangleupdate);
            if sense==1
                x_new =  round(x_previous+x_step); y_new =  round(y_previous+y_step); z_new =  round(z_previous+z_step);
            else
                x_new =  round(x_previous-x_step); y_new =  round(y_previous-y_step); z_new =  round(z_previous-z_step);
            end
        end

        % Not out of bounds
        x_new = min([x_new domain_size(1)]); y_new = min([y_new domain_size(2)]); z_new = min([z_new domain_size(3)]);
        x_new = max([x_new 1]); y_new = max([y_new 1]); z_new = max([z_new 1]);

        % Generate candidate tube
        n_step = max(abs([x_new, y_new, z_new] - [x_previous, y_previous, z_previous]))+1;
        idx = round((x_new-x_previous)/(n_step-1) * [0:1:n_step-1] + x_previous);
        idy = round((y_new-y_previous)/(n_step-1) * [0:1:n_step-1] + y_previous);
        idz = round((z_new-z_previous)/(n_step-1) * [0:1:n_step-1] + z_previous);
        idQ = sub2ind(domain_size,idx,idy,idz);
        tube_tmp = zeros(domain_size);
        tube_tmp(idQ)=1;
        % Add volume
        if radius>0
            dmap = bwdist(tube_tmp,p.shape);
            tube_tmp(dmap<=radius)=1;
        end

        % Check overlapping
        overlapping_tube = 1; overlaping_mask = 1;
        if ~p.tube_can_overlap_with_tube
            overlapping_tube = tube_tmp + tube_phase;
        end
        if ~tube_can_overlap_with_mask
            overlaping_mask = tube_tmp + ~p.maskused;
        end
        if max(max(max(overlapping_tube)))==1 && max(max(max(overlaping_mask)))==1 % No overlapping with other tube and with mask
            % We can add the tube
            points = {[x_previous y_previous z_previous]; [x_new y_new z_new] };
            length_tube = sqrt(x_step^2 + y_step^2 + z_step^2);

            % Update stoping condition
            k_tube = k_tube + 1;
            if p.seeds_onlyatedges && p.force_seedequidistribution
                seed_xyz = seed_xyz+seed_candidate;
            end
            volume_fraction = sum(sum(sum(tube_phase~=0)))/nvoxel;
            % Update results
            tube_skeleton(idQ) = k_tube;
            idx_tube = find(tube_tmp);
            tube_id(idx_tube) = k_tube;
            tube_phase(idx_tube) = 1;
            tubes(k_tube).radius = radius;
            tubes(k_tube).points = points;
            % Update seed location (for verification only)
            used_seeds(x_seed,y_seed,z_seed)=1;

            % Update seeds
            seeds = seeds .* (~tube_phase);
            idseeds = find(seeds);
            if sav.save_continuously
                function_save_tif(uint8(tube_phase),[sav.Savefolder 'Tubelabel_run_' num2str(p.k_run) '_tube' num2str(k_tube) '.tif']);
                if k_tube<=255
                    function_save_tif(uint8(tube_id),[sav.Savefolder 'TubeId_run_' num2str(p.k_run) '_tube' num2str(k_tube) '.tif']);
                    function_save_tif(uint8(tube_skeleton),[sav.Savefolder 'Tubeskeleton_run_ ' num2str(p.k_run) '_tube' num2str(k_tube) '.tif']);
                else
                    function_save_tif(uint16(tube_id),[sav.Savefolder 'TubeId_run_' num2str(p.k_run) '_tube' num2str(k_tube) '.tif']);
                    function_save_tif(uint16(tube_skeleton),[sav.Savefolder 'Tubeskeleton_run_ ' num2str(p.k_run) '_tube' num2str(k_tube) '.tif']);
                end
                save([sav.Savefolder 'TubeInfo_run_' num2str(p.k_run) '_tube' num2str(k_tube) '.mat'],'tubes','-mat');
            end
            fprintf('- Tube #%i/%i generated with %i control points and length %1.3f. Volume fraction %1.5f/%1.5f\n',k_tube,stopcond.n_tube,2,length_tube(end),volume_fraction,stopcond.target_volumefraction)

            if check_isotropy
                if mod(k_tube,10)==0
                    k_iso = k_iso+1;
                    
                    % [tau_pore, aniso_pore] = calculate_isotropy(~tube_phase);
                    % if k_iso==1
                    %     x_iso = k_tube;
                    %     tau_pore_fig = [tau_pore(1) tau_pore(2) tau_pore(3)];
                    %     aniso_pore_fig = aniso_pore;
                    %     Fig_pore = figure;
                    %     ax_pore=axes('Parent',Fig_pore);
                    % else
                    %     x_iso = [x_iso;k_tube];
                    %     tau_pore_fig = [tau_pore_fig; [tau_pore(1) tau_pore(2) tau_pore(3)]];
                    %     aniso_pore_fig = [aniso_pore_fig; aniso_pore];
                    %     cla(ax_pore);
                    % end
                    % hold(ax_pore,"on");
                    % plot(x_iso,tau_pore_fig(:,1),'Parent',ax_pore,'Color','b');
                    % plot(x_iso,tau_pore_fig(:,2),'Parent',ax_pore,'Color','r');
                    % plot(x_iso,tau_pore_fig(:,3),'Parent',ax_pore,'Color','g');
                    % plot(x_iso,aniso_pore_fig,'Parent',ax_pore,'Color','k');
                    % hold(ax_pore,"off");
                    % tau_pore_fig

                    [tau_tube, aniso_tube] = calculate_isotropy(tube_phase);
                    if k_iso==1
                        x_iso = k_tube;
                        tau_tube_fig = [tau_tube(1) tau_tube(2) tau_tube(3)];
                        aniso_tube_fig = aniso_tube;
                        Fig_tube = figure;
                        ax_tube=axes('Parent',Fig_tube);
                    else
                        x_iso = [x_iso;k_tube];
                        tau_tube_fig = [tau_tube_fig; [tau_tube(1) tau_tube(2) tau_tube(3)]];
                        aniso_tube_fig = [aniso_tube_fig; aniso_tube];
                        cla(ax_tube);
                    end
                    hold(ax_tube,"on");
                    plot(x_iso,tau_tube_fig(:,1),'Parent',ax_tube,'Color','b');
                    plot(x_iso,tau_tube_fig(:,2),'Parent',ax_tube,'Color','r');
                    plot(x_iso,tau_tube_fig(:,3),'Parent',ax_tube,'Color','g');
                    plot(x_iso,aniso_tube_fig,'Parent',ax_tube,'Color','k');
                    hold(ax_tube,"off");
                    tau_tube_fig
                    pause(0.1)
                end
            end

            % Reset timer
            tStart_tube = tic;
            continue % Move to next tube

        else
            time_currenttube = toc(tStart_tube);
            continue % Pick another seed and another orientation
        end

    else % Curved tube, more complicated code

        % Initialize
        points = cell([1,1e4]);
        length_tube = zeros(1,1e4);
        npoint = 1;
        points(npoint) = {[x_previous y_previous z_previous]};

        % Reset try attempt
        n_attempt = 0;
        n_trygoback = 0;
        goback = false;
        within_goback_attempt = false;
        n_point_start = [];

        reached_edge = false;
        while true % Exit only through break statment

            if goback
                if n_trygoback==1
                    n_point_start = npoint-1;
                    n_point_end = npoint+1;
                    within_goback_attempt = true;
                end
                points(npoint) =  {[]};
                npoint = npoint-1;
                tmp = cell2mat(points(npoint));
                x_previous = tmp(1);
                y_previous = tmp(2);
                z_previous = tmp(3);
                n_attempt = 0;
                goback = false;
            end
            tmp_points = points;

            if npoint>1 % New orientation
                if strcmp(p.change_orientation_choice,'Angle')
                    azimuth_change = pickfromdistriubtion(p.maxdeviation_azimuthchange,p.equiprobability_azimuthchange,p.mean_azimuthchange,p.cdf_azimuthchange,p.x_azimuthchange);
                    elevation_change = pickfromdistriubtion(p.maxdeviation_elevationchange,p.equiprobability_elevationchange,p.mean_elevationchange,p.cdf_elevationchange,p.x_elevationchange);
                else
                    Ux_change = pickfromdistriubtion(p.maxdeviation_Uxchange,p.equiprobability_Uxchange,p.mean_Uxchange,p.cdf_Uxchange,p.x_Uxchange);
                    Uy_change = pickfromdistriubtion(p.maxdeviation_Uychange,p.equiprobability_Uychange,p.mean_Uychange,p.cdf_Uychange,p.x_Uychange);
                    Uz_change = pickfromdistriubtion(p.maxdeviation_Uzchange,p.equiprobability_Uzchange,p.mean_Uzchange,p.cdf_Uzchange,p.x_Uzchange);
                    %r = rand(3,1);
                    %u = r/norm(r);
                    %[azimuth_change,elevation_change,~] = cart2sph(u(1),u(2),u(3));
                    [azimuth_change,elevation_change,~] = cart2sph(Ux_change,Uy_change,Uz_change);
                    azimuth_change =rad2deg(azimuth_change);
                    elevation_change =rad2deg(elevation_change);                    
                end
                if p.relaxation_orientation
                    azimuth_change = azimuth_change + n_attempt/p.relaxation_rate * azimuth_change;
                    elevation_change = elevation_change + n_attempt/p.relaxation_rate * elevation_change;
                end
                azimuth = azimuth + azimuth_change;
                elevation = elevation + elevation_change;
            end

            % Distance to next control point
            distancebetweenangleupdate = pickfromdistriubtion(p.maxdeviation_distancebetweenangleupdate,p.equiprobability_distancebetweenangleupdate,p.mean_distancebetweenangleupdate,p.cdf_distancebetweenangleupdate,p.x_distancebetweenangleupdate);

            [x_step,y_step,z_step] = sph2cart2(deg2rad(azimuth),deg2rad(elevation),distancebetweenangleupdate);
            distance = sqrt(x_step^2 + y_step^2 + z_step^2);
            x_new =  round(x_previous+x_step);
            y_new =  round(y_previous+y_step);
            z_new =  round(z_previous+z_step);

            if x_new<1 || y_new<1 || z_new<1 || x_new>domain_size(1) || y_new>domain_size(2) || z_new>domain_size(3) % Out of bound
                if npoint==1 % We need another seed
                    if p.verbose
                        fprintf('Control point out of bound (1 point) -> discard tube\n')
                    end
                    time_currenttube = toc(tStart_tube);
                    break
                else % put point at the eddge
                    slopes = [(x_new-x_previous) (y_new-y_previous) (z_new-z_previous)]/distancebetweenangleupdate;
                    backward = [(1-x_previous) (1-y_previous) (1-z_previous)]./slopes;
                    forward = [(domain_size(1)-x_previous) (domain_size(2)-y_previous) (domain_size(3)-z_previous)]./slopes;
                    backward(backward<0)=0;
                    forward(forward<0)=0;
                    t = backward+forward;
                    idx = t==min(t);
                    distancebetweenangleupdate = t(idx)+1; 
                    distancebetweenangleupdate = distancebetweenangleupdate(1); % In case of equality
                    [x_step,y_step,z_step] = sph2cart2(deg2rad(azimuth),deg2rad(elevation),distancebetweenangleupdate);
                    x_new =  round(x_previous+x_step); y_new =  round(y_previous+y_step); z_new =  round(z_previous+z_step);

                    % Not out of bounds
                    x_new = min([x_new domain_size(1)]); y_new = min([y_new domain_size(2)]); z_new = min([z_new domain_size(3)]);
                    x_new = max([x_new 1]); y_new = max([y_new 1]); z_new = max([z_new 1]);

                    % Verification
                    if x_new~=1 && y_new~=1 && z_new~=1 && x_new~=domain_size(1) && y_new~=domain_size(2) && z_new~=domain_size(3)
                        warning('Control point has not been put at the edges!')
                        keyboard
                    end

                end
            end

            if x_new==1 || y_new==1 || z_new==1 || x_new==domain_size(1) || y_new==domain_size(2) || z_new==domain_size(3) % point at the eddge
                reached_edge = true; % If last candidate point is valid, and distance requirement is ok, then keep the tube
            end


            % Check if going in the wrong direction
            if p.check_backwardprogression && npoint>p.backwardprogression_maxnpoint
                s = cell2mat(points(1));
                p0 =  cell2mat(points(npoint-p.backwardprogression_maxnpoint+1));
                p1 =  cell2mat(points(npoint));
                dp0 = sum((p0-s).^2).^0.5;
                dp1 = sum((p1-s).^2).^0.5;
                if dp1<dp0
                    if p.verbose
                        fprintf('Wrong direction -> discard tube\n')
                    end
                    break % Try with another seed
                end
            end

            % Grow tube
            tmp_points(npoint+1) = {[x_new y_new z_new]};
            tmp_points(npoint+2:end)=[];
            point2point_distance = sqrt(x_step^2 + y_step^2 + z_step^2);
            current_length = length_tube(npoint) + point2point_distance;

            % check contiguity
            Q =[0 0 0;2 2 2];
            iterprecision = 0;
            while max(abs(diff( Q(:,1) ))) >1 || max(abs(diff( Q(:,2) ))) >1 || max(abs(diff( Q(:,3) ))) >1
                iterprecision = iterprecision+1;
                % % CALL THIRD PARTY CODE
                Q = custom_hobbysplines(tmp_points,'debug',false,'cycle',Bezier_pars.cycle,'tension',Bezier_pars.tension,'Nq',round(current_length*Qprecision*iterprecision));
                % % END OF THIRD PARTY CODE
                Q = round(Q);
            end


            if min(Q(:,1)) >= 1 && min(Q(:,2)) >= 1 && min(Q(:,3)) >= 1 && max(Q(:,1)) <= domain_size(1) && max(Q(:,2)) <= domain_size(2) && max(Q(:,3)) <= domain_size(3) && ~anynan(Q) % Tube is within bounds
                idQ = sub2ind(domain_size,Q(:,1),Q(:,2),Q(:,3));
                idQ(isnan(idQ))=[];
                % Create tube
                tube_tmp = zeros(domain_size);
                tube_tmp(idQ) = 1;

                if radius>0
                    dmap = bwdist(tube_tmp,p.shape);
                    tube_tmp(dmap<=radius)=1;
                end

                % Check overlapping
                overlapping_tube = 1; overlaping_mask = 1;
                if ~p.tube_can_overlap_with_tube
                    overlapping_tube = tube_tmp + tube_phase;
                end
                if ~tube_can_overlap_with_mask
                    overlaping_mask = tube_tmp + ~p.maskused;
                end
                if max(max(max(overlapping_tube)))==1 && max(max(max(overlaping_mask)))==1 % No overlapping with other tube and with mask
                    % Add point
                    npoint = npoint+1;
                    if p.verbose
                        fprintf('   Current number of point %i\n',npoint);
                    end
                    %npoint
                    points(npoint) = tmp_points(npoint);
                    x_previous = x_new;
                    y_previous = y_new;
                    z_previous = z_new;
                    n_attempt = 0; % reset
                    if within_goback_attempt && npoint==n_point_end
                        n_trygoback=0; % Reset
                        within_goback_attempt = false;
                        n_point_start = [];
                    end
                    length_tube(npoint) = length_tube(npoint-1)+point2point_distance;


                    % Check tube
                    if (p.finitelength && length_tube(npoint)>=tubelength) || (~p.finitelength && reached_edge && length_tube(npoint)>=p.min_length) % Keep tube

                        % Update stoping condition
                        k_tube = k_tube + 1;
                        volume_fraction = sum(sum(sum(tube_phase~=0)))/nvoxel;
                        if p.seeds_onlyatedges && p.force_seedequidistribution
                            seed_xyz = seed_xyz+seed_candidate;
                        end
                        % Update results
                        tube_skeleton(idQ) = k_tube;
                        idx_tube = find(tube_tmp);
                        tube_id(idx_tube) = k_tube;
                        tube_phase(idx_tube) = 1;
                        tubes(k_tube).radius = radius;
                        points(npoint+1:end)=[];
                        tubes(k_tube).points = points;
                        % Update seed location (for verification only)
                        used_seeds(x_seed,y_seed,z_seed)=1;
                        % Update seeds
                        seeds = seeds .* (~tube_phase);
                        idseeds = find(seeds);
                        if sav.save_continuously
                            function_save_tif(uint8(tube_phase),[sav.Savefolder 'Tubelabel_run_' num2str(p.k_run) '_tube' num2str(k_tube) '.tif']);
                            if k_tube<=255
                                function_save_tif(uint8(tube_id),[sav.Savefolder 'TubeId_run_' num2str(p.k_run) '_tube' num2str(k_tube) '.tif']);
                                function_save_tif(uint8(tube_skeleton),[sav.Savefolder 'Tubeskeleton_run_ ' num2str(p.k_run) '_tube' num2str(k_tube) '.tif']);
                            else
                                function_save_tif(uint16(tube_id),[sav.Savefolder 'TubeId_run_' num2str(p.k_run) '_tube' num2str(k_tube) '.tif']);
                                function_save_tif(uint16(tube_skeleton),[sav.Savefolder 'Tubeskeleton_run_ ' num2str(p.k_run) '_tube' num2str(k_tube) '.tif']);
                            end
                            save([sav.Savefolder 'TubeInfo_run_' num2str(p.k_run) '_tube' num2str(k_tube) '.mat'],'tubes','-mat');
                        end
                        fprintf('- Tube #%i/%i generated with %i control points and length %1.3f. Volume fraction %1.5f/%1.5f\n',k_tube,stopcond.n_tube,npoint,length_tube(npoint),volume_fraction,stopcond.target_volumefraction)
                        % Reset timer
                        tStart_tube = tic;
                        break % Move to next tube
                    end
                    if reached_edge && ( (p.finitelength && length_tube(npoint)<tubelength) || (~p.finitelength && length_tube(npoint)<p.min_length) ) % Discard tube
                        if p.verbose
                            fprintf('Reached edge but too short -> discard tube\n')
                        end
                        time_currenttube = toc(tStart_tube);
                        break
                    end

                else % Overlapping
                    n_attempt = n_attempt+1;
                    if n_attempt <= p.n_attempt_max
                        if p.verbose
                            fprintf('      Overlapping -> Try different orientation\n');
                        end
                        continue % Try with another orientation from same starting point
                    else
                        n_trygoback = n_trygoback + 1;
                        if n_trygoback <= p.n_trygoback_max && npoint>2 && (~within_goback_attempt || (within_goback_attempt && (npoint-1)==n_point_start))
                            goback = true; % Try with another orientation, but starting from the previous point
                            if p.verbose
                                fprintf('      Overlapping. Dead end case -> Move back to previous point\n');
                            end
                            continue
                        else
                            if p.verbose
                                fprintf('      Overlapping. Dead end case -> try other seed\n');
                            end
                            break % Try with another seed
                        end
                    end
                end

            else % While tube control points are in bounds, tube Bezier curve is out of bound
                if npoint==1 % We need another seed
                    if p.verbose
                        fprintf('Tube out of bound (1 point) -> discard tube\n')
                    end
                    break
                else
                    n_attempt = n_attempt+1;
                    if n_attempt <= p.n_attempt_max
                        if p.verbose
                            fprintf('      Tube out of bound -> try different orientation\n');
                        end
                        continue % Try with another orientation from same starting point
                    else
                        n_trygoback = n_trygoback + 1;
                        if n_trygoback <= p.n_trygoback_max && npoint>2 && (~within_goback_attempt || (within_goback_attempt && npoint-1==n_point_start))
                            goback = true; % Try with another orientation, but starting from the previous point
                            if p.verbose
                                fprintf('      Tube out of bound. Dead end case -> Move back to previous point\n');
                            end
                            continue
                        else
                            break % Try with another seed
                        end
                    end
                end
            end

        end

    end

end

if sum(cropp.cut)>0
    l1 = cropp.cut(1);
    l2 = cropp.cut(2);
    l3 = cropp.cut(3);
    tube_phase = tube_phase(l1:end-l1+1, l2:end-l2+1, l3:end-l3+1);
    tube_id = tube_id(l1:end-l1+1, l2:end-l2+1, l3:end-l3+1);
    tube_skeleton = tube_skeleton(l1:end-l1+1, l2:end-l2+1, l3:end-l3+1);
end

if visualize_seeds
    idx = find(used_seeds);
    [IX,IY,IZ]=ind2sub(domain_size,idx);
    sum(IX==1)
    sum(IX==domain_size(1))
    sum(IY==1)
    sum(IY==domain_size(2))
    sum(IZ==1)
    sum(IZ==domain_size(3))
    Microstructure_basic_visualization_interface(used_seeds)
end

%% LOCAL FUNCTION
    function [val] = pickfromdistriubtion(maxdeviation,equiprobability,mean_val,cdf_,x_val)
        if maxdeviation==0
            val = mean_val;
        else
            if equiprobability
                val = rand*(2*maxdeviation) + (mean_val-maxdeviation);
            else
                val = interp1(cdf_,x_val,rand);
            end
        end
    end

    function [tau, aniso] = calculate_isotropy(BW)
        tau = zeros(1,3);
        Tau_factor_result = TauFactor('InLine',1,0,0,BW,[0 0 0;0 0 0;1 0 0],[1 1 1]);
        if Tau_factor_result.Tau_W1.Tau == Inf % No percolation path. 65535 is default value for inf tortuosity factor
            Tau_factor_result.Tau_W1.Tau = NaN;
        end
        tau(1) = Tau_factor_result.Tau_W1.Tau;

        Tau_factor_result = TauFactor('InLine',1,0,0,BW,[0 0 0;0 0 0;0 1 0],[1 1 1]);
        if Tau_factor_result.Tau_W2.Tau == Inf
            Tau_factor_result.Tau_W2.Tau = NaN;
        end
        tau(2) = Tau_factor_result.Tau_W2.Tau;

        Tau_factor_result = TauFactor('InLine',1,0,0,BW,[0 0 0;0 0 0;0 0 1],[1 1 1]);
        if Tau_factor_result.Tau_W3.Tau == Inf
            Tau_factor_result.Tau_W3.Tau = NaN;
        end
        tau(3) = Tau_factor_result.Tau_W3.Tau;
             
        aniso = nanmax(tau) - nanmin(tau);
    end

end