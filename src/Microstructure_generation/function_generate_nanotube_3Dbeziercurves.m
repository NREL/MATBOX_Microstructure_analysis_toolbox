function [tube_phase, tube_id, tube_skeleton, tubes, outcome] = function_generate_nanotube_3Dbeziercurves(domain_size,Bezier_pars,sav,stopcond,cropp,p)

%% THIRD-PARTY LICENSE
% I am using a slighly modified function "hobbysplines" to create smooth 3D bezier curves from:

% Will Robertson (2023). Smooth 3D bezier curves with implicit control points
% https://www.mathworks.com/matlabcentral/fileexchange/42302-smooth-3d-bezier-curves-with-implicit-control-points
% MATLAB Central File Exchange. Retrieved March 29, 2023.
% License is available within hobbysplines.m

% Everything else is NREL code covered with the MATBOX licence

%% HARD CODED PARAMETERS
check_contiguity = true;
check_contiguity_eachtube = false; % tube non continguity is rare, so only check at the end to save time
Nq = 5000; % Set high to avoid non contiguity

%% TUBE ORIENTATION, ANGULAR DEVIATION, RADIUS, and UPDATE DISTANCE
% Cumulative and distribution functions
[x_initialazimuth, ~, cdf_initialazimuth] = distributions(p.equiprobability_initialazimuth, p.mean_initialazimuth, p.sigma_initialazimuth, p.maxdeviation_initialazimuth);
[x_initialelevation, ~, cdf_initialelevation] = distributions(p.equiprobability_initialelevation, p.mean_initialelevation, p.sigma_initialelevation, p.maxdeviation_initialelevation);
[x_azimuthchange, ~, cdf_azimuthchange] = distributions(p.equiprobability_azimuthchange, p.mean_azimuthchange, p.sigma_azimuthchange, p.maxdeviation_azimuthchange);
[x_elevationchange, ~, cdf_elevationchange] = distributions(p.equiprobability_elevationchange, p.mean_elevationchange, p.sigma_elevationchange, p.maxdeviation_elevationchange);
[x_distancebetweenangleupdate, ~, cdf_distancebetweenangleupdate] = distributions(p.equiprobability_distancebetweenangleupdate, p.mean_distancebetweenangleupdate, p.sigma_distancebetweenangleupdate, p.maxdeviation_distancebetweenangleupdate);
[x_radius, ~, cdf_radius] = distributions(p.equiprobability_radius, p.mean_radius, p.sigma_radius, p.maxdeviation_radius);
if p.finitelength
    [x_length, ~, cdf_length] = distributions(p.equiprobability_length, p.mean_length, p.sigma_length, p.maxdeviation_length);
end

%% GENERATION
% Initialization
tube_phase = zeros(domain_size);
tube_id = zeros(domain_size);
tube_skeleton = zeros(domain_size);
if isempty(p.maskused)
    p.maskused = zeros(domain_size);
end

%% SET TUBE SEEDS
% =1 possible starting location
% =0 forbidden starting location

% Initialization
edges = tube_phase+1;

if p.seeds_onlyatedges % Locate seeds at domain's edge
    edges(2:end-1,2:end-1,2:end-1)=0;
    if sum(cropp.distfromedge)>0
        l1 = round(domain_size(1)*cropp.distfromedge(1));
        l2 = round(domain_size(2)*cropp.distfromedge(2));
        l3 = round(domain_size(3)*cropp.distfromedge(3));

        edges(:,1:l2,1:l3)=0; edges(:,1:l2,end-l3+1:end)=0;
        edges(:,end-l2+1:end,1:l3)=0; edges(:,end-l2+1:end,end-l3+1:end)=0;

        edges(1:l1,:,1:l3)=0; edges(1:l1,:,end-l3+1:end)=0;
        edges(end-l1+1:end,:,1:l3)=0; edges(end-l1+1:end,:,end-l3+1:end)=0;

        edges(1:l1,1:l2,:)=0; edges(1:l1,end-l2+1:end,:)=0;
        edges(end-l1+1:end,1:l2,:)=0; edges(end-l1+1:end,end-l2+1:end,:)=0;
    end
end

if ~isempty(p.maskused) % Starting location can not be within mask
    edges = edges .* (~p.maskused);
end

idseeds = find(edges); % Index of all initial starting locations

%% ALGORITHM
outcome = 'Success (#tube)'; % By default
tStart_total = tic;   
tStart_tube = tic;   

nvoxel = numel(tube_phase);
k_tube = 0;
volume_fraction = 0;
tau = [];
gotoedge = false;
while k_tube < p.n_tube
    % Stoping condition: volume fraction
    if volume_fraction>=p.target_volumefraction
        outcome = 'Success (vf)'; % By default
        return;
    end
    % Stoping condition: time
    t_total = toc(tStart_total);
    t_tube = toc(tStart_tube);
    if t_total > stopcond.totaltime
        outcome = 'Exceed total time';
        return;
    end
    if t_tube > stopcond.pertube
        outcome = 'Exceed tube time';
        return;
    end

    % Select radius within a distribution
    if p.maxdeviation_radius==0
        radius = p.mean_radius;
    else
        if p.equiprobability_radius
            radius = rand*(2*p.maxdeviation_radius) + (p.mean_radius-p.maxdeviation_radius);
        else
            radius = interp1(cdf_radius,x_radius,rand);
        end
    end

    % Select length within a distribution
    if p.finitelength
        if p.maxdeviation_length==0
            tubelength = p.mean_length;
        else
            if p.equiprobability_length
                tubelength = rand*(2*p.maxdeviation_length) + (p.mean_length-p.maxdeviation_length);
            else
                tubelength = interp1(cdf_length,x_length,rand);
            end
        end
    end

    % Select initial orientation
    if p.maxdeviation_initialazimuth==0
        azimuth = p.mean_initialazimuth;
    else
        if p.equiprobability_initialazimuth
            azimuth = rand*(2*p.maxdeviation_initialazimuth) + (p.mean_initialazimuth-p.maxdeviation_initialazimuth);
        else
            azimuth = interp1(cdf_initialazimuth,x_initialazimuth,rand);
        end
    end
    if p.maxdeviation_initialelevation==0
        elevation = p.mean_initialelevation;
    else
        if p.equiprobability_initialelevation
            elevation = rand*(2*p.maxdeviation_initialelevation) + (p.mean_initialelevation-p.maxdeviation_initialelevation);
        else
            elevation = interp1(cdf_initialelevation,x_initialelevation,rand);
        end
    end

    % Select seed
    idseed = max([1 round(numel(idseeds)*rand)]);
    [x_previous, y_previous, z_previous] = ind2sub(domain_size,idseeds(idseed));

    % Initialize
    npoint = 1;
    points = cell([1,1e6]);
    points(npoint) = {[x_previous y_previous z_previous]};
    length_tube = zeros(1,1e6);

    n_attempt = 0; % reset
    goback = false;
    n_trygoback = 0;
    within_goback_attempt = false;
    n_point_start = [];
    first_outofbound = true;
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
            if p.maxdeviation_azimuthchange==0
                azimuth_change = p.mean_azimuthchange;
            else
                if p.equiprobability_azimuthchange
                    azimuth_change = rand*(2*p.maxdeviation_azimuthchange) + (p.mean_azimuthchange-p.maxdeviation_azimuthchange);
                else
                    azimuth_change = interp1(cdf_azimuthchange,x_azimuthchange,rand);
                end
            end
            if p.maxdeviation_elevationchange==0
                elevation_change = p.mean_elevationchange;
            else
                if p.equiprobability_elevationchange
                    elevation_change = rand*(2*p.maxdeviation_elevationchange) + (p.mean_elevationchange-p.maxdeviation_elevationchange);
                else
                    elevation_change = interp1(cdf_elevationchange,x_elevationchange,rand);
                end
            end

azimuth_change
elevation_change

            if p.relaxation_orientation
                azimuth_change = azimuth_change + n_attempt/p.relaxation_rate * azimuth_change;
                elevation_change = elevation_change + n_attempt/p.relaxation_rate * elevation_change;
            end

            azimuth = azimuth + azimuth_change;
            elevation = elevation + elevation_change;
            
        end
        % Distance to next control point
        if p.maxdeviation_distancebetweenangleupdate==0
            distancebetweenangleupdate = p.mean_distancebetweenangleupdate;
        else
            if p.equiprobability_distancebetweenangleupdate
                distancebetweenangleupdate = rand*(2*p.maxdeviation_distancebetweenangleupdate) + (p.mean_distancebetweenangleupdate-p.maxdeviation_distancebetweenangleupdate);
            else
                distancebetweenangleupdate = interp1(cdf_distancebetweenangleupdate,x_distancebetweenangleupdate,rand);
            end
        end

        if p.maxdeviation_azimuthchange==0 && p.maxdeviation_elevationchange==0 % Let's go directly to end point
            if p.finitelength
                distancebetweenangleupdate = 1.001*tubelength;
            else
                distancebetweenangleupdate = 1000; % Force scaling step length
            end
            gotoedge = true;
        end

        [xstep,ystep,zstep] = sph2cart(deg2rad(azimuth),deg2rad(elevation),distancebetweenangleupdate);
        distance = sqrt(xstep^2 + ystep^2 + zstep^2);
        xnew =  round(x_previous+xstep);
        ynew =  round(y_previous+ystep);
        znew =  round(z_previous+zstep);

        [xnew ynew znew]
        keyboard

        if xnew<1 || ynew<1 || znew<1 || xnew>domain_size(1) || ynew>domain_size(2) || znew>domain_size(3) % Out of bounds
            if npoint==1 && ~gotoedge % We need another seed
                break
            else
                if (p.finitelength && length_tube(npoint)<tubelength) || (~p.finitelength && length_tube(npoint)<p.min_length) % Discard, not long enough. Try with another seed
                    break
                else
                    if first_outofbound % Put point at the edge
                        backward = [(1-x_previous)/xstep (1-y_previous)/ystep (1-z_previous)/zstep];
                        cond1 = backward>0;
                        cond2 = backward<1;
                        id1 = find(cond1.*cond2);
                        k1  = [];
                        if ~isempty(id1)
                            k1 = min(backward(id1));
                        end

                        forward = [(domain_size(1)-x_previous)/xstep (domain_size(2)-y_previous)/ystep (domain_size(3)-z_previous)/zstep];
                        cond1 = forward>0;
                        cond2 = forward<1;
                        id2 = find(cond1.*cond2);
                        k2 =[];
                        if ~isempty(id2)
                            k2 = min(forward(id2));
                        end    
                        scale_step = min([k1 k2]);

                        [xstep,ystep,zstep] = sph2cart(deg2rad(azimuth),deg2rad(elevation),scale_step*distancebetweenangleupdate);
                        distance = sqrt(xstep^2 + ystep^2 + zstep^2);
                        xnew =  round(x_previous+xstep);
                        ynew =  round(y_previous+ystep);
                        znew =  round(z_previous+zstep);

                        x_new = min([xnew domain_size(1)]); y_new = min([ynew domain_size(2)]); z_new = min([znew domain_size(3)]);
                        x_new = max([xnew 1]); y_new = max([ynew 1]); z_new = max([znew 1]);
                        
                        first_outofbound = false;

                    else % Keep tube
                        k_tube = k_tube + 1;

                        % tube_id = tube_id + tube_tmp*k_tube; % Sometimes, there is a tiny overlapp between tubes
                        idx_tube = find(tube_tmp);
                        tube_id(idx_tube) = k_tube;

                        tube_skeleton(idQ) = k_tube;

                        tube_phase = tube_phase + tube_tmp;
                        tube_phase(tube_phase~=0)=1; % Sometimes, there is a tiny overlapp between tubes

                        tubes(k_tube).radius = radius;
                        points(npoint+1:end)=[];
                        tubes(k_tube).points = points;

                        volume_fraction = sum(sum(sum(tube_phase~=0)))/nvoxel;

                        % Update seed
                        edges = edges .* (~tube_phase);
                        idseeds = find(edges);

                        if sav.save_continuously
                            % function_save_tif(uint8(tube_phase),[sav.Savefolder 'Tubelabel_run_' num2str(p.k_run) '.tif']);
                            % if k_tube<=255
                            %    function_save_tif(uint8(tube_id),[sav.Savefolder 'TubeId_run_' num2str(p.k_run) '.tif']);
                            %    function_save_tif(uint8(tube_skeleton),[sav.Savefolder 'Tubeskeleton_run_ ' num2str(p.k_run) '.tif']);
                            % else
                            %     function_save_tif(uint16(tube_id),[sav.Savefolder 'TubeId_run_' num2str(p.k_run) '.tif']);
                            %     function_save_tif(uint16(tube_skeleton),[sav.Savefolder 'Tubeskeleton_run_ ' num2str(p.k_run) '.tif']);
                            % end
                            % save([sav.Savefolder 'TubeInfo_run_' num2str(p.k_run) '.mat'],'tubes','-mat');
                            % Above may have issue with write permission

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
                        fprintf('- Tube #%i/%i generated with %i control points and length %1.3f. Volume fraction %1.5f/%1.5f\n',k_tube,p.n_tube,npoint,length_tube(npoint),volume_fraction,p.target_volumefraction)
                        tStart_tube = tic; % rest timer for new tube

                        if check_contiguity
                            if check_contiguity_eachtube || (~check_contiguity_eachtube && (volume_fraction>=p.target_volumefraction || k_tube>=p.n_tube))

                                disp 'Checking tube contiguity';
                                ids = unique(tube_id);
                                ids(ids==0)=[]; % remove background
                                find_nocontiguous = false;
                                for kid=1:1:length(ids)
                                    BW = zeros(size(tube_phase));
                                    idx = find(tube_id==ids(kid));
                                    BW(idx)=1;
                                    L = bwlabeln(BW,26);
                                    if length(unique(L))>2 % Non-contiguous
                                        tube_phase(idx)=0;
                                        tube_id(idx)=0;
                                        tube_skeleton(tube_skeleton==ids(kid)) = 0;
                                        find_nocontiguous = true;
                                        fprintf('- Tube #%i is not continguous and has been removed\n',ids(kid));
                                    end
                                end
                                if find_nocontiguous
                                    disp 'Renumerotating...';
                                    % Recalculate volume fraction
                                    volume_fraction = sum(sum(sum(tube_phase~=0)))/nvoxel;
                                    % Renumeroate
                                    old_ids = unique(tube_id);
                                    old_ids(old_ids==0)=[];
                                    k_tube = length(old_ids);
                                    new_ids = 1:1:k_tube;
                                    tmp = tube_id;
                                    for kid=1:1:k_tube
                                        tube_id(tmp == old_ids(kid)) = new_ids(kid);
                                    end
                                    clear tmp;
                                    tmp = tube_skeleton;
                                    for kid=1:1:k_tube
                                        tube_skeleton(tmp == old_ids(kid)) = new_ids(kid);
                                    end
                                    clear tmp;

                                    tmp = tubes;
                                    clear tubes;
                                    for kid=1:1:k_tube
                                        radius = tmp(old_ids(kid)).radius;
                                        points = tmp(old_ids(kid)).points;
                                        tubes(kid).radius = radius;
                                        tubes(kid).points = points;
                                    end

                                else
                                    disp 'All tubes are contiguous';
                                end

                            end
                        end


                        break

                    end
                end

            end
        end


        if isempty(xnew)
            warning('step is empty')
            keyboard
        end

        % Check if within bounds
        if xnew>=1 && ynew>=1 && znew>=1 && xnew<=domain_size(1) && ynew<=domain_size(2) && znew<=domain_size(3)
            % Check if goind wrong direction
            if p.check_backwardprogression && npoint>p.backwardprogression_maxnpoint
                s = cell2mat(points(1));
                p0 =  cell2mat(points(npoint-p.backwardprogression_maxnpoint+1));
                p1 =  cell2mat(points(npoint));
                dp0 = sum((p0-s).^2).^0.5;
                dp1 = sum((p1-s).^2).^0.5;
                if dp1<dp0
                    %disp 'wrong direction'
                    break % Try with another seed
                end
            end
                        
            tmp_points(npoint+1) = {[xnew ynew znew]};
            tmp_points(npoint+2:end)=[];

            % % CALL THIRD PARTY CODE
            Q = custom_hobbysplines(tmp_points,'debug',false,'cycle',Bezier_pars.cycle,'tension',Bezier_pars.tension,'Nq',Nq);
            % % END OF THIRD PARTY CODE

            Q = round(Q);
            if min(Q(:,1)) >= 1 && min(Q(:,2)) >= 1 && min(Q(:,3)) >= 1 && max(Q(:,1)) <= domain_size(1) && max(Q(:,2)) <= domain_size(3) && max(Q(:,3)) <= domain_size(3)
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
                overlapping1 = tube_tmp + tube_phase;
                overlapping2 = tube_tmp + p.maskused;
                %overlapping = tube_tmp + tube_phase + p.maskused;
                % if max(max(max(overlapping)))==1 % No overlapping with other tube and with mask. Last tube point can be within mask!
                if max(max(max(overlapping1)))==1 && max(max(max(overlapping2)))==1 % No overlapping with other tube and with mask
                    % Add point
                    npoint = npoint+1
                    %npoint
                    points(npoint) = tmp_points(npoint);
                    x_previous = xnew;
                    y_previous = ynew;
                    z_previous = znew;
                    n_attempt = 0; % reset
                    if within_goback_attempt && npoint==n_point_end
                        n_trygoback=0; % Reset
                        within_goback_attempt = false;
                        n_point_start = [];
                    end
                    length_tube(npoint) = length_tube(npoint-1)+distance;

                    if p.finitelength && length_tube(npoint)>=tubelength
                        % Keep
                        k_tube = k_tube + 1;

                        % tube_id = tube_id + tube_tmp*k_tube; % Sometimes, there is a tiny overlapp between tubes
                        idx_tube = find(tube_tmp);
                        tube_id(idx_tube) = k_tube;

                        tube_skeleton(idQ) = k_tube;

                        tube_phase = tube_phase + tube_tmp;
                        tube_phase(tube_phase~=0)=1; % Sometimes, there is a tiny overlapp between tubes

                        tubes(k_tube).radius = radius;
                        points(npoint+1:end)=[];
                        tubes(k_tube).points = points;

                        volume_fraction = sum(sum(sum(tube_phase~=0)))/nvoxel;

                        % Update seed
                        edges = edges .* (~tube_phase);
                        idseeds = find(edges);

                        if sav.save_continuously
                            % function_save_tif(uint8(tube_phase),[sav.Savefolder 'Tubelabel_run_' num2str(p.k_run) '.tif']);
                            % if k_tube<=255
                            %    function_save_tif(uint8(tube_id),[sav.Savefolder 'TubeId_run_' num2str(p.k_run) '.tif']);
                            %    function_save_tif(uint8(tube_skeleton),[sav.Savefolder 'Tubeskeleton_run_ ' num2str(p.k_run) '.tif']);
                            % else
                            %     function_save_tif(uint16(tube_id),[sav.Savefolder 'TubeId_run_' num2str(p.k_run) '.tif']);
                            %     function_save_tif(uint16(tube_skeleton),[sav.Savefolder 'Tubeskeleton_run_ ' num2str(p.k_run) '.tif']);
                            % end
                            % save([sav.Savefolder 'TubeInfo_run_' num2str(p.k_run) '.mat'],'tubes','-mat');
                            % Above may have issue with write permission

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
                        fprintf('- Tube #%i/%i generated with %i control points and length %1.3f. Volume fraction %1.5f/%1.5f\n',k_tube,p.n_tube,npoint,length_tube(npoint),volume_fraction,p.target_volumefraction)
                        tStart_tube = tic; % rest timer for new tube

                        if check_contiguity
                            if check_contiguity_eachtube || (~check_contiguity_eachtube && (volume_fraction>=p.target_volumefraction || k_tube>=p.n_tube))

                                disp 'Checking tube contiguity';
                                ids = unique(tube_id);
                                ids(ids==0)=[]; % remove background
                                find_nocontiguous = false;
                                for kid=1:1:length(ids)
                                    BW = zeros(size(tube_phase));
                                    idx = find(tube_id==ids(kid));
                                    BW(idx)=1;
                                    L = bwlabeln(BW,26);
                                    if length(unique(L))>2 % Non-contiguous
                                        tube_phase(idx)=0;
                                        tube_id(idx)=0;
                                        tube_skeleton(tube_skeleton==ids(kid)) = 0;
                                        find_nocontiguous = true;
                                        fprintf('- Tube #%i is not continguous and has been removed\n',ids(kid));
                                    end
                                end
                                if find_nocontiguous
                                    disp 'Renumerotating...';
                                    % Recalculate volume fraction
                                    volume_fraction = sum(sum(sum(tube_phase~=0)))/nvoxel;
                                    % Renumeroate
                                    old_ids = unique(tube_id);
                                    old_ids(old_ids==0)=[];
                                    k_tube = length(old_ids);
                                    new_ids = 1:1:k_tube;
                                    tmp = tube_id;
                                    for kid=1:1:k_tube
                                        tube_id(tmp == old_ids(kid)) = new_ids(kid);
                                    end
                                    clear tmp;
                                    tmp = tube_skeleton;
                                    for kid=1:1:k_tube
                                        tube_skeleton(tmp == old_ids(kid)) = new_ids(kid);
                                    end
                                    clear tmp;

                                    tmp = tubes;
                                    clear tubes;
                                    for kid=1:1:k_tube
                                        radius = tmp(old_ids(kid)).radius;
                                        points = tmp(old_ids(kid)).points;
                                        tubes(kid).radius = radius;
                                        tubes(kid).points = points;
                                    end

                                else
                                    disp 'All tubes are contiguous';
                                end

                            end
                        end




                        break

                    end


                else
                    n_attempt = n_attempt+1;
                    if n_attempt <= p.n_attempt_max
                        continue % Try with another orientation from same starting point
                    else
                        n_trygoback = n_trygoback + 1;

                        % n_trygoback
                        % within_goback_attempt
                        % npoint-1
                        % if ~isempty(n_point_start)
                        %     n_point_start
                        % end
                        

                        if n_trygoback <= p.n_trygoback_max && npoint>2 && (~within_goback_attempt || (within_goback_attempt && (npoint-1)==n_point_start))
                            goback = true; % Try with another orientation, but starting from the previous point
                            %[npoint n_trygoback]
                            continue
                        else
                            %disp 'new seed';
                            %keyboard
                            break % Try with another seed
                        end
                    end
                end

            else
                if npoint==1 % We need another seed
                    break
                else
                    n_attempt = n_attempt+1;
                    if n_attempt <= p.n_attempt_max
                        continue % Try with another orientation from same starting point
                    else
                        n_trygoback = n_trygoback + 1;
                        if n_trygoback <= p.n_trygoback_max && npoint>2 && (~within_goback_attempt || (within_goback_attempt && npoint-1==n_point_start))
                            goback = true; % Try with another orientation, but starting from the previous point
                            %[npoint n_trygoback]
                            continue
                        else
                            break % Try with another seed
                        end
                    end
                end
            end

        else % Out of bounds
           warning('This should not be reached!')
        end

        % At the edges
        if npoint>1 && (xnew==1 || ynew==1 || znew==1 || xnew==domain_size(1) || ynew==domain_size(2) || znew==domain_size(3))
            if (p.finitelength && length_tube(npoint)<tubelength) || (~p.finitelength && length_tube(npoint)<p.min_length) % Discard, not long enough. Try with another seed
                break
            else
                % Keep
                k_tube = k_tube + 1;
                % tube_id = tube_id + tube_tmp*k_tube; % Sometimes, there is a tiny overlapp between tubes
                idx_tube = find(tube_tmp);
                tube_id(idx_tube) = k_tube;
                tube_skeleton(idQ) = k_tube;
                tube_phase = tube_phase + tube_tmp;
                tube_phase(tube_phase~=0)=1; % Sometimes, there is a tiny overlapp between tubes

                tubes(k_tube).radius = radius;
                points(npoint+1:end)=[];
                tubes(k_tube).points = points;

                volume_fraction = sum(sum(sum(tube_phase~=0)))/nvoxel;

                % Update seed
                edges = edges .* (~tube_phase);
                idseeds = find(edges);

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

                fprintf('- Tube #%i/%i generated with %i control points and length %1.3f. Volume fraction %1.5f/%1.5f\n',k_tube,p.n_tube,npoint,length_tube(npoint),volume_fraction,p.target_volumefraction)
                tStart_tube = tic; % rest timer for new tube

                if check_contiguity
                    if check_contiguity_eachtube || (~check_contiguity_eachtube && (volume_fraction>=p.target_volumefraction || k_tube>=p.n_tube))
                        
                        disp 'Checking tube contiguity';
                        ids = unique(tube_id);
                        ids(ids==0)=[]; % remove background
                        find_nocontiguous = false;
                        for kid=1:1:length(ids)
                            BW = zeros(size(tube_phase));
                            idx = find(tube_id==ids(kid));
                            BW(idx)=1;
                            L = bwlabeln(BW,26);
                            if length(unique(L))>2 % Non-contiguous
                                tube_phase(idx)=0;
                                tube_id(idx)=0;
                                tube_skeleton(tube_skeleton==ids(kid)) = 0;
                                find_nocontiguous = true;
                                fprintf('- Tube #%i is not continguous and has been removed\n',ids(kid));
                            end
                        end 
                        if find_nocontiguous
                            disp 'Renumerotating...';
                            % Recalculate volume fraction
                            volume_fraction = sum(sum(sum(tube_phase~=0)))/nvoxel;
                            % Renumeroate
                            old_ids = unique(tube_id);
                            old_ids(old_ids==0)=[];
                            k_tube = length(old_ids);
                            new_ids = 1:1:k_tube;
                            tmp = tube_id;
                            for kid=1:1:k_tube
                                tube_id(tmp == old_ids(kid)) = new_ids(kid);
                            end
                            clear tmp;
                            tmp = tube_skeleton;
                            for kid=1:1:k_tube
                                tube_skeleton(tmp == old_ids(kid)) = new_ids(kid);
                            end
                            clear tmp;

                            tmp = tubes;
                            clear tubes;
                            for kid=1:1:k_tube
                                radius = tmp(old_ids(kid)).radius;
                                points = tmp(old_ids(kid)).points;
                                tubes(kid).radius = radius;
                                tubes(kid).points = points;
                            end                            

                        else
                            disp 'All tubes are contiguous';
                        end

                    end
                end


                break
            end

        end

    end

    % % Solid tortuosity check
    % [ntau,~] = size(tau);
    % if k_tube > ntau
    % 
    % 
    %     Tau_factor_result = TauFactor('InLine',1,0,0,tube_phase,[0 0 0;0 0 0;1 0 0],[1 1 1]);
    %     if Tau_factor_result.Tau_W1.Tau == Inf % No percolation path. 65535 is default value for inf tortuosity factor
    %         Tau_factor_result.Tau_W1.Tau = NaN;
    %     end
    %     Tau_factor_result = TauFactor('InLine',1,0,0,tube_phase,[0 0 0;0 0 0;0 1 0],[1 1 1]);
    %     if Tau_factor_result.Tau_W2.Tau == Inf
    %         Tau_factor_result.Tau_W2.Tau = NaN;
    %     end
    %     Tau_factor_result = TauFactor('InLine',1,0,0,tube_phase,[0 0 0;0 0 0;0 0 1],[1 1 1]);
    %     if Tau_factor_result.Tau_W3.Tau == Inf
    %         Tau_factor_result.Tau_W3.Tau = NaN;
    %     end
    % 
    %     Tau_factor_result.Tau_W1.Tau
    %     tau(k_tube,1) = Tau_factor_result.Tau_W1.Tau;
    %     tau(k_tube,2) = Tau_factor_result.Tau_W2.Tau;
    %     tau(k_tube,3) = Tau_factor_result.Tau_W3.Tau;
    %     aniso(k_tube) = nanmax(tau(k_tube,:)) - nanmin(tau(k_tube,:));
    % 
    %     tau
    % 
    %     if k_tube==1
    %         Fig = figure;
    %         ax_=axes('Parent',Fig);
    %     else
    %         cla(ax_);
    %     end
    %     plot([1:1:k_tube],aniso(k_tube))




    %end
                    

end

if sum(cropp.cut)>0
    l1 = cropp.cut(1);
    l2 = cropp.cut(2);
    l3 = cropp.cut(3);
    tube_phase = tube_phase(l1:end-l1+1, l2:end-l2+1, l3:end-l3+1);
    tube_id = tube_id(l1:end-l1+1, l2:end-l2+1, l3:end-l3+1);
    tube_skeleton = tube_skeleton(l1:end-l1+1, l2:end-l2+1, l3:end-l3+1);
end


%% LOCAL FUNCTION
    function [x, pdf_, cdf_] = distributions(equiprobability, mean_, sigma_, maxdeviation)
        x = []; pdf_ = []; cdf_ =[];

        % Univariate Gaussian Distribution
        % Probability density function of a gaussian distribution
        pdf_GD = @(x,mu,sigma) 1./(2*pi*sigma.^2).^(0.5).*exp(-(x-mu).^2 ./ (2*sigma.^2));

        if ~equiprobability && maxdeviation~=0
            x=linspace(mean_-maxdeviation, mean_+maxdeviation, 100);
            pdf_=pdf_GD(x,mean_,sigma_);
            cdf_ = pdf2cdf(x,pdf_);
        end

        function c = pdf2cdf(x,pdf)
            n=length(x);
            c=zeros(1,n);
            for k=2:1:n
                c(k) = trapz(x(1:k),pdf(1:k));
            end
            if c(end)<1
                c(end)=1;
            end
        end
    end  

end
