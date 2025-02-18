function [All_subdomain,GROUP_SUBDOMAIN, REP2ndsize, Representativity_2ndlength_str, Representativity_analysis_size2] = RVE_Subdomains(pRVE, sz, voxel_size)
% Determine the bounds of subdomains from a domain.

All_subdomain = [];
GROUP_SUBDOMAIN = [];
REP2ndsize = [];
Representativity_2ndlength_str = []; % empty, or 'square root', or 'length';
Representativity_analysis_size2 = [];

if length(sz)==2 || sz(3)==1
    dimension = 2;
    sz(3) = 1;
else
    dimension = 3;
end

divisions = pRVE.divisions;

% Initialize
sub_domain_id=0;
sub_domain_group=0;

if strcmp(pRVE.type,'A') % 'Independant subvolumes of same size + keep initial aspect ratio (A)'
    initial_volume = prod(sz);
    initial_volume_um = initial_volume*voxel_size^dimension;
    if pRVE.donotcut_smaller
        max_division = floor((initial_volume_um/(pRVE.donotcut_smaller_than^dimension))^(1/dimension));
        divisions(divisions>max_division) = [];
    end

    if ~isempty(divisions)
        for k=1:1:length(divisions)
            sub_domain_group=sub_domain_group+1;
            number_subdomain_group=divisions(k)^dimension;
            GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
            if dimension == 3
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (sz(1)/divisions(k) * sz(2)/divisions(k) * sz(3)/divisions(k))^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1)/divisions(k) sz(2)/divisions(k) sz(3)/divisions(k)] ./ (sz(3)/divisions(k));
                zk = divisions(k);
            else
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (sz(1)/divisions(k) * sz(2)/divisions(k) )^(1/2);
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1)/divisions(k) sz(2)/divisions(k) ] ./ (sz(2)/divisions(k));
                zk = 1;
            end
            tmp_id = [];
            for x_=1:1:divisions(k)
                for y_=1:1:divisions(k)
                    for z_=1:1:zk

                        x0 = floor((sz(1)/divisions(k)*(x_-1))+1);
                        x1 = floor(sz(1)/divisions(k)*x_);
                        y0 = floor((sz(2)/divisions(k)*(y_-1))+1);
                        y1 = floor(sz(2)/divisions(k)*y_);
                        if dimension == 3
                            z0 = floor((sz(3)/divisions(k)*(z_-1))+1);
                            z1 = floor(sz(3)/divisions(k)*z_);
                        else
                            z0=1;
                            z1=1;
                        end

                        sub_domain_id = sub_domain_id+1;
                        All_subdomain(sub_domain_id,1) = sub_domain_id;
                        All_subdomain(sub_domain_id,2) = sub_domain_group;
                        All_subdomain(sub_domain_id,3) = number_subdomain_group;
                        All_subdomain(sub_domain_id,4) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length;
                        All_subdomain(sub_domain_id,5) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length;
                        All_subdomain(sub_domain_id,6) = GROUP_SUBDOMAIN.id(sub_domain_group).length;                        
                        All_subdomain(sub_domain_id,7) = x0;
                        All_subdomain(sub_domain_id,8) = x1;
                        All_subdomain(sub_domain_id,9) = y0;
                        All_subdomain(sub_domain_id,10) = y1;
                        All_subdomain(sub_domain_id,11) = z0;
                        All_subdomain(sub_domain_id,12) = z1;
                        All_subdomain(sub_domain_id,13) = divisions(k);

                        tmp_id = [tmp_id sub_domain_id];
                    end
                end
            end
            GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
        end
    end    
end


if  strcmp(pRVE.type,'B') % 'Independant subvolumes of same size + user-defined aspect ratio (B)'
    
    if dimension==2
        AR_desired = pRVE.Aspectratio_RVE / pRVE.Aspectratio_RVE(2);
        % Search for the largest subvolume that has the desired AR
        l1_n = 1; l1(1) = sz(1)/l1_n;
        l2(1) = AR_desired(2)*l1(1)/AR_desired(1); l2_n = sz(2)/l2(1);
        vol = [0 0];
        if l2_n>=1
            vol(1)=prod([l1(1) l2(1)]);
        end

        l2_n = 1; l2(2) = sz(2)/l2_n;
        l1(2) = AR_desired(1)*l2(2)/AR_desired(2); l1_n = sz(1)/l1(2);
        if l1_n>=1
            vol(2)=prod([l1(2) l2(2)]);
        end
    else
        AR_desired = pRVE.Aspectratio_RVE / pRVE.Aspectratio_RVE(3);
        % Search for the largest subvolume that has the desired AR
        l1_n = 1; l1(1) = sz(1)/l1_n;
        l2(1) = AR_desired(2)*l1(1)/AR_desired(1); l2_n = sz(2)/l2(1);
        l3(1) = AR_desired(3)*l1(1)/AR_desired(1); l3_n = sz(3)/l3(1);
        vol = [0 0 0];
        if l2_n>=1 && l3_n>=1
            vol(1)=prod([l1(1) l2(1) l3(1)]);
        end

        l2_n = 1; l2(2) = sz(2)/l2_n;
        l1(2) = AR_desired(1)*l2(2)/AR_desired(2); l1_n = sz(1)/l1(2);
        l3(2) = AR_desired(3)*l2(2)/AR_desired(2); l3_n = sz(3)/l3(2);
        if l1_n>=1 && l3_n>=1
            vol(2)=prod([l1(2) l2(2) l3(2)]);
        end

        l3_n = 1; l3(3) = sz(3)/l3_n;
        l1(3) = AR_desired(1)*l3(3)/AR_desired(3); l1_n = sz(1)/l1(3);
        l2(3) = AR_desired(2)*l3(3)/AR_desired(3); l2_n = sz(2)/l2(3);
        if l1_n>=1 && l2_n>=1
            vol(3)=prod([l1(3) l2(3) l3(3)]);
        end
    end

    idx=find(vol==max(vol)); % Case that ensure the largest volume
    idx=idx(1); % In case of equality
    if dimension==2
        Max_subdomain_size = [l1(idx) l2(idx)];
    else
        Max_subdomain_size = [l1(idx) l2(idx) l3(idx)];
    end

    initial_volume = prod(Max_subdomain_size);
    initial_volume_um = initial_volume*voxel_size^dimension;
    if pRVE.donotcut_smaller
        max_division = floor((initial_volume_um/(pRVE.donotcut_smaller_than^dimension))^(1/dimension));
        divisions(divisions>max_division) = [];
    end    

    if ~isempty(divisions)
        for k=0:1:length(divisions)
            sub_domain_group=sub_domain_group+1;
            if k==0
                current_subdomain_size = Max_subdomain_size;
            else
                current_subdomain_size = Max_subdomain_size/divisions(k);
            end
            div_ = floor(sz./current_subdomain_size);
            current_domain_size = div_.*current_subdomain_size;
            number_subdomain_group=prod(div_);

            GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
            if dimension == 2
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (current_domain_size(1)/div_(1) * current_domain_size(2)/div_(2))^(1/2);
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [current_domain_size(1)/div_(1) current_domain_size(2)/div_(2)] ./ (current_domain_size(2)/div_(2));
                div_(3) = 1;
            else
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (current_domain_size(1)/div_(1) * current_domain_size(2)/div_(2) * current_domain_size(3)/div_(3))^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [current_domain_size(1)/div_(1) current_domain_size(2)/div_(2) current_domain_size(3)/div_(3)] ./ (current_domain_size(3)/div_(3));
            end
            tmp_id = [];
            for x_=1:1:div_(1)
                for y_=1:1:div_(2)
                    for z_=1:1:div_(3)

                        x0 = floor((current_domain_size(1)/div_(1)*(x_-1))+1);
                        x1 = floor(current_domain_size(1)/div_(1)*x_);
                        y0 = floor((current_domain_size(2)/div_(2)*(y_-1))+1);
                        y1 = floor(current_domain_size(2)/div_(2)*y_);
                        if dimension == 2
                            z0 = 1;
                            z1 = 1;
                        else
                            z0 = floor((current_domain_size(3)/div_(3)*(z_-1))+1);
                            z1 = floor(current_domain_size(3)/div_(3)*z_);
                        end

                        sub_domain_id = sub_domain_id+1;
                        actual_size = ((x1-x0+1)*(y1-y0+1)*(z1-z0+1))^(1/3);
                        All_subdomain(sub_domain_id,1) = sub_domain_id;
                        All_subdomain(sub_domain_id,2) = sub_domain_group;
                        All_subdomain(sub_domain_id,3) = number_subdomain_group;
                        All_subdomain(sub_domain_id,4) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length;
                        All_subdomain(sub_domain_id,5) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length;
                        All_subdomain(sub_domain_id,6) = GROUP_SUBDOMAIN.id(sub_domain_group).length;
                        All_subdomain(sub_domain_id,7) = x0;
                        All_subdomain(sub_domain_id,8) = x1;
                        All_subdomain(sub_domain_id,9) = y0;
                        All_subdomain(sub_domain_id,10) = y1;
                        All_subdomain(sub_domain_id,11) = z0;
                        All_subdomain(sub_domain_id,12) = z1;
                        All_subdomain(sub_domain_id,13) = 0;

                        tmp_id = [tmp_id sub_domain_id];
                    end
                end
            end
            GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
        end
    end
end


if strcmp(pRVE.type,'C') % 'Independant subvolumes of same size + one constant length (C)'
    
    if strcmp(pRVE.Constantdirection,'Direction 1')
        initial_section = sz(2)*sz(3);
    elseif strcmp(pRVE.Constantdirection,'Direction 2')
        initial_section = sz(1)*sz(3);
    elseif strcmp(pRVE.Constantdirection,'Direction 3')
        initial_section = sz(1)*sz(2);
    end    
    initial_section_um = initial_section*voxel_size^(dimension-1);
    if pRVE.donotcut_smaller
        max_division = floor( (initial_section_um/(pRVE.donotcut_smaller_than^(dimension-1)))^(1/(dimension-1)) );
        divisions(divisions>max_division) = [];
    end

    if dimension == 3
        Representativity_2ndlength_str = 'square root';
    else
        Representativity_2ndlength_str = 'length';
    end

    if ~isempty(divisions)
        for k=1:1:length(divisions)
            sub_domain_group=sub_domain_group+1; % Increment group id
            number_subdomain_group=divisions(k)^(dimension-1); % Number of subdomain
            if strcmp(pRVE.Constantdirection,'Direction 1')
                REP2ndsize = (sz(2) * sz(3))^(1/(dimension-1));
                if dimension == 3
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (sz(2)/divisions(k) * sz(3)/divisions(k) * sz(1))^(1/3);
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (sz(2)/divisions(k) * sz(3)/divisions(k))^(1/2);
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1) sz(2)/divisions(k) sz(3)/divisions(k)] ./ (sz(3)/divisions(k));
                    zk = divisions(k);                    
                else
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (sz(1) * sz(2)/divisions(k))^(1/2);
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = sz(2)/divisions(k);
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1) sz(2)/divisions(k)] ./ (sz(2)/divisions(k));
                    zk = 1;
                end
                tmp_id = [];
                for y_=1:1:divisions(k)
                    for z_=1:1:zk
                        x0 = 1;
                        x1 = sz(1);
                        y0 = floor((sz(2)/divisions(k)*(y_-1))+1);
                        y1 = floor(sz(2)/divisions(k)*y_);
                        if dimension == 3
                            z0 = floor((sz(3)/divisions(k)*(z_-1))+1);
                            z1 = floor(sz(3)/divisions(k)*z_);
                        else
                            z0 = 1;
                            z1 = 1;
                        end
                        sub_domain_id = sub_domain_id+1;

                        All_subdomain(sub_domain_id,1) = sub_domain_id;
                        All_subdomain(sub_domain_id,2) = sub_domain_group;
                        All_subdomain(sub_domain_id,3) = number_subdomain_group;
                        All_subdomain(sub_domain_id,4) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length;
                        All_subdomain(sub_domain_id,5) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length;
                        All_subdomain(sub_domain_id,6) = GROUP_SUBDOMAIN.id(sub_domain_group).length;
                        All_subdomain(sub_domain_id,7) = x0;
                        All_subdomain(sub_domain_id,8) = x1;
                        All_subdomain(sub_domain_id,9) = y0;
                        All_subdomain(sub_domain_id,10) = y1;
                        All_subdomain(sub_domain_id,11) = z0;
                        All_subdomain(sub_domain_id,12) = z1;
                        All_subdomain(sub_domain_id,13) = divisions(k);

                        tmp_id = [tmp_id sub_domain_id];
                    end
                end

            elseif strcmp(pRVE.Constantdirection,'Direction 2')
                REP2ndsize = (sz(1) * sz(3))^(1/(dimension-1));
                if dimension == 3
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (sz(1)/divisions(k) * sz(3)/divisions(k) * sz(2))^(1/3);
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (sz(1)/divisions(k) * sz(3)/divisions(k))^(1/2);
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1)/divisions(k) sz(2) sz(3)/divisions(k)] ./ (sz(3)/divisions(k));
                    zk = divisions(k);
                else
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (sz(1)/divisions(k) * sz(2))^(1/2);
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = sz(1)/divisions(k);
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1)/divisions(k) sz(2)/divisions(k)] ./ (sz(2));
                    zk = 1;
                end
                tmp_id = [];
                for x_=1:1:divisions(k)
                    for z_=1:1:zk
                        x0 = floor((sz(1)/divisions(k)*(x_-1))+1);
                        x1 = floor(sz(1)/divisions(k)*x_);
                        y0 = 1;
                        y1 = sz(2);
                        if dimension == 3
                            z0 = floor((sz(3)/divisions(k)*(z_-1))+1);
                            z1 = floor(sz(3)/divisions(k)*z_);
                        else
                            z0 = 1;
                            z1 = 1;
                        end
                        sub_domain_id = sub_domain_id+1;

                        All_subdomain(sub_domain_id,1) = sub_domain_id;
                        All_subdomain(sub_domain_id,2) = sub_domain_group;
                        All_subdomain(sub_domain_id,3) = number_subdomain_group;
                        All_subdomain(sub_domain_id,4) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length;
                        All_subdomain(sub_domain_id,5) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length;
                        All_subdomain(sub_domain_id,6) = GROUP_SUBDOMAIN.id(sub_domain_group).length;
                        All_subdomain(sub_domain_id,7) = x0;
                        All_subdomain(sub_domain_id,8) = x1;
                        All_subdomain(sub_domain_id,9) = y0;
                        All_subdomain(sub_domain_id,10) = y1;
                        All_subdomain(sub_domain_id,11) = z0;
                        All_subdomain(sub_domain_id,12) = z1;
                        All_subdomain(sub_domain_id,13) = divisions(k);

                        tmp_id = [tmp_id sub_domain_id];
                    end
                end

            elseif strcmp(pRVE.Constantdirection,'Direction 3')
                if dimension == 3
                    REP2ndsize = (sz(1) * sz(2))^(1/2);
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (sz(1)/divisions(k) * sz(2)/divisions(k) * sz(3))^(1/3);
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (sz(1)/divisions(k) * sz(2)/divisions(k))^(1/2);
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1)/divisions(k) sz(2)/divisions(k) sz(3)] ./ sz(3);
                    tmp_id = [];
                    for x_=1:1:divisions(k)
                        for y_=1:1:divisions(k)
                            x0 = floor((sz(1)/divisions(k)*(x_-1))+1);
                            x1 = floor(sz(1)/divisions(k)*x_);
                            y0 = floor((sz(2)/divisions(k)*(y_-1))+1);
                            y1 = floor(sz(2)/divisions(k)*y_);
                            z0 = 1;
                            z1 = sz(3);

                            sub_domain_id = sub_domain_id+1;
                            %actual_size = ((x1-x0+1)*(y1-y0+1)*(z1))^(1/3);
                            %actual_section = ((x1-x0+1)*(y1-y0+1))^(1/2);

                            All_subdomain(sub_domain_id,1) = sub_domain_id;
                            All_subdomain(sub_domain_id,2) = sub_domain_group;
                            All_subdomain(sub_domain_id,3) = number_subdomain_group;
                            All_subdomain(sub_domain_id,4) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length;
                            All_subdomain(sub_domain_id,5) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length;
                            All_subdomain(sub_domain_id,6) = GROUP_SUBDOMAIN.id(sub_domain_group).length;
                            All_subdomain(sub_domain_id,7) = x0;
                            All_subdomain(sub_domain_id,8) = x1;
                            All_subdomain(sub_domain_id,9) = y0;
                            All_subdomain(sub_domain_id,10) = y1;
                            All_subdomain(sub_domain_id,11) = z0;
                            All_subdomain(sub_domain_id,12) = z1;
                            All_subdomain(sub_domain_id,13) = divisions(k);

                            tmp_id = [tmp_id sub_domain_id];
                        end
                    end
                else
                    warning('Impossible RVE analysis: domain is 2D but constant direction is set as Direction 3');
                end
            end
            GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
        end
    end

    if pRVE.subs2

        if pRVE.donotcut_smaller
            if (initial_section_um/2)^(1/(dimension-1)) > pRVE.donotcut_smaller_than

                sub_domain_group=sub_domain_group+1; % Increment group id
                number_subdomain_group=2; % Number of subdomain
                tmp_id = [];

                if dimension == 3
                    add2 = true;
                    if strcmp(pRVE.Constantdirection,'Direction 1')
                        if sz(2)>=sz(3)
                            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (sz(2)/2 * sz(3) * sz(1))^(1/3);
                            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (sz(2)/2 * sz(3))^(1/2);
                            GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                            GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1) sz(2)/2 sz(3)] ./ sz(3);
                            yy_start(1)=1; yy_end(1)=floor(sz(2)/2);
                            yy_start(2)=yy_end(1)+1; yy_end(2)=sz(2);
                            zz_start(1)=1; zz_end(1)=sz(3);
                            zz_start(2)=1; zz_end(2)=sz(3);
                        else
                            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (sz(2) * sz(3)/2 * sz(1))^(1/3);
                            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (sz(2) * sz(3)/2)^(1/2);
                            GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                            GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1) sz(2) sz(3)/2] ./ (sz(3)/2);
                            yy_start(1)=1; yy_end(1)=sz(2);
                            yy_start(2)=1; yy_end(2)=sz(2);
                            zz_start(1)=1; zz_end(1)=floor(sz(3)/2);
                            zz_start(2)=zz_end(1)+1; zz_end(2)=sz(3);
                        end
                        xx_start(1) = 1; xx_end(1) = sz(1);
                        xx_start(2) = 1; xx_end(2) = sz(1);

                    elseif strcmp(pRVE.Constantdirection,'Direction 2')
                        if sz(1)>=sz(3)
                            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (sz(1)/2 * sz(3) * sz(2))^(1/3);
                            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (sz(1)/2 * sz(3))^(1/2);
                            GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                            GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1)/2 sz(2) sz(3)] ./ sz(3);
                            xx_start(1)=1; xx_end(1)=floor(sz(1)/2);
                            xx_start(2)=xx_end(1)+1; xx_end(2)=sz(1);
                            zz_start(1)=1; zz_end(1)=sz(3);
                            zz_start(2)=1; zz_end(2)=sz(3);
                        else
                            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (sz(1) * sz(3)/2 * sz(2))^(1/3);
                            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (sz(1) * sz(3)/2)^(1/2);
                            GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                            GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1) sz(2) sz(3)/2] ./ (sz(3)/2);
                            xx_start(1)=1; xx_end(1)=sz(1);
                            xx_start(2)=1; xx_end(2)=sz(1);
                            zz_start(1)=1; zz_end(1)=floor(sz(3)/2);
                            zz_start(2)=zz_end(1)+1; zz_end(2)=sz(3);
                        end
                        yy_start(1) = 1; yy_end(1) = sz(2);
                        yy_start(2) = 1; yy_end(2) = sz(2);

                    elseif strcmp(pRVE.Constantdirection,'Direction 3')
                        if sz(1)>=sz(2)
                            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (sz(1)/2 * sz(2) * sz(3))^(1/3);
                            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (sz(1)/2 * sz(2))^(1/2);
                            GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                            GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1)/2 sz(2) sz(3)] ./ sz(3);
                            xx_start(1)=1; xx_end(1)=floor(sz(1)/2);
                            xx_start(2)=xx_end(1)+1; xx_end(2)=sz(1);
                            yy_start(1)=1; yy_end(1)=sz(2);
                            yy_start(2)=1; yy_end(2)=sz(2);
                        else
                            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (sz(1) * sz(2)/2 * sz(3))^(1/3);
                            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (sz(1) * sz(2)/2)^(1/2);
                            GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                            GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1) sz(2)/2 sz(3)] ./ sz(3);
                            xx_start(1)=1; xx_end(1)=sz(1);
                            xx_start(2)=1; xx_end(2)=sz(1);
                            yy_start(1)=1; yy_end(1)=floor(sz(2)/2);
                            yy_start(2)=yy_end(1)+1; yy_end(2)=sz(2);
                        end
                        zz_start(1) = 1; zz_end(1) = sz(3);
                        zz_start(2) = 1; zz_end(2) = sz(3);
                    end
                else
                    if sum(divisions==2)==0
                        add2 = true;
                        if strcmp(pRVE.Constantdirection,'Direction 1')
                            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = 0;
                            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (sz(1) * sz(2)/2)^(1/2);
                            GROUP_SUBDOMAIN.id(sub_domain_group).length = sz(2)/2;
                            GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1) sz(2)/2] ./ (sz(2)/2);
                            xx_start(1) = 1; xx_end(1) = sz(1);
                            xx_start(2) = 1; xx_end(2) = sz(1);
                            yy_start(1)=1; yy_end(1)=floor(sz(2)/2);
                            yy_start(2)=yy_end(1)+1; yy_end(2)=sz(2);
                            zz_start(1)=1; zz_end(1)=1;
                            zz_start(2)=1; zz_end(2)=1;

                        elseif strcmp(pRVE.Constantdirection,'Direction 2')
                            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = 0;
                            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (sz(1)/2 * sz(2))^(1/2);
                            GROUP_SUBDOMAIN.id(sub_domain_group).length = sz(1)/2;
                            GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1)/2 sz(2)] ./ (sz(2));
                            xx_start(1)=1; xx_end(1)=floor(sz(1)/2);
                            xx_start(2)=xx_end(1)+1; yy_end(2)=sz(1);
                            yy_start(1) = 1; yy_end(1) = sz(2);
                            yy_start(2) = 1; yy_end(2) = sz(2);
                            zz_start(1)=1; zz_end(1)=1;
                            zz_start(2)=1; zz_end(2)=1;
                        elseif strcmp(pRVE.Constantdirection,'Direction 3')
                            warning('Impossible RVE analysis: domain is 2D but constant direction is set as Direction 3');
                        end
                    else
                        add2 = false; % No need, it has been already done
                    end
                end

                if add2
                    for sub2=1:1:2
                        sub_domain_id = sub_domain_id+1;
                        All_subdomain(sub_domain_id,1) = sub_domain_id;
                        All_subdomain(sub_domain_id,2) = sub_domain_group;
                        All_subdomain(sub_domain_id,3) = number_subdomain_group;
                        All_subdomain(sub_domain_id,4) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length;
                        All_subdomain(sub_domain_id,5) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length;
                        All_subdomain(sub_domain_id,6) = GROUP_SUBDOMAIN.id(sub_domain_group).length;
                        All_subdomain(sub_domain_id,7) = xx_start(sub2);
                        All_subdomain(sub_domain_id,8) = xx_end(sub2);
                        All_subdomain(sub_domain_id,9) = yy_start(sub2);
                        All_subdomain(sub_domain_id,10) = yy_end(sub2);
                        All_subdomain(sub_domain_id,11) = zz_start(sub2);
                        All_subdomain(sub_domain_id,12) = zz_end(sub2);
                        All_subdomain(sub_domain_id,13) = 0;
                        tmp_id = [tmp_id sub_domain_id];
                    end
                    GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
                end

            end
        end

    end
end

if strcmp(pRVE.type,'D') % 'Independant subvolumes of same size + two constant lengths (D)'
    
    if strcmp(pRVE.Constantdirection,'Directions 1 and 2')
        initial_length = sz(3);
    elseif strcmp(pRVE.Constantdirection,'Direction 1 and 3')
        initial_length = sz(2);
    elseif strcmp(pRVE.Constantdirection,'Direction 2 and 3')
        initial_length = sz(1);
    end      
    initial_length_um = initial_length*voxel_size;
    if pRVE.donotcut_smaller
        max_division = floor( initial_length_um/pRVE.donotcut_smaller_than );
        divisions(divisions>max_division) = [];
    end
    Representativity_2ndlength_str = 'length';

    if ~isempty(divisions)
        for k=1:1:length(divisions)
            sub_domain_group=sub_domain_group+1; % Increment group id
            number_subdomain_group=divisions(k)^1; % Number of subdomain
            if strcmp(pRVE.Constantdirection,'Directions 1 and 2')
                REP2ndsize = sz(3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = ( sz(1) * sz(2) * sz(3)/divisions(k) )^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).length = sz(3)/divisions(k);
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1) sz(2) sz(3)/divisions(k)] ./ (sz(3)/divisions(k));
                tmp_id = [];
                for z_=1:1:divisions(k)
                    x0 = 1;
                    x1 = sz(1);
                    y0 = 1;
                    y1 = sz(2);
                    z0 = floor((sz(3)/divisions(k)*(z_-1))+1);
                    z1 = floor(sz(3)/divisions(k)*z_);

                    sub_domain_id = sub_domain_id+1;

                    All_subdomain(sub_domain_id,1) = sub_domain_id;
                    All_subdomain(sub_domain_id,2) = sub_domain_group;
                    All_subdomain(sub_domain_id,3) = number_subdomain_group;
                    All_subdomain(sub_domain_id,4) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length;
                    All_subdomain(sub_domain_id,5) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length;
                    All_subdomain(sub_domain_id,6) = GROUP_SUBDOMAIN.id(sub_domain_group).length;
                    All_subdomain(sub_domain_id,7) = x0;
                    All_subdomain(sub_domain_id,8) = x1;
                    All_subdomain(sub_domain_id,9) = y0;
                    All_subdomain(sub_domain_id,10) = y1;
                    All_subdomain(sub_domain_id,11) = z0;
                    All_subdomain(sub_domain_id,12) = z1;
                    All_subdomain(sub_domain_id,13) = divisions(k);

                    tmp_id = [tmp_id sub_domain_id];
                end

            elseif strcmp(pRVE.Constantdirection,'Directions 1 and 3')
                REP2ndsize = sz(2);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = ( sz(1) * sz(2)/divisions(k) * sz(3) )^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).length = sz(2)/divisions(k);
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1) sz(2)/divisions(k) sz(3)] ./ sz(3);
                tmp_id = [];
                for y_=1:1:divisions(k)
                    x0 = 1;
                    x1 = sz(1);
                    y0 = floor((sz(2)/divisions(k)*(y_-1))+1);
                    y1 = floor(sz(2)/divisions(k)*y_);
                    z0 = 1;
                    z1 = sz(3);

                    sub_domain_id = sub_domain_id+1;

                    All_subdomain(sub_domain_id,1) = sub_domain_id;
                    All_subdomain(sub_domain_id,2) = sub_domain_group;
                    All_subdomain(sub_domain_id,3) = number_subdomain_group;
                    All_subdomain(sub_domain_id,4) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length;
                    All_subdomain(sub_domain_id,5) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length;
                    All_subdomain(sub_domain_id,6) = GROUP_SUBDOMAIN.id(sub_domain_group).length;
                    All_subdomain(sub_domain_id,7) = x0;
                    All_subdomain(sub_domain_id,8) = x1;
                    All_subdomain(sub_domain_id,9) = y0;
                    All_subdomain(sub_domain_id,10) = y1;
                    All_subdomain(sub_domain_id,11) = z0;
                    All_subdomain(sub_domain_id,12) = z1;
                    All_subdomain(sub_domain_id,13) = divisions(k);

                    tmp_id = [tmp_id sub_domain_id];
                end

            elseif strcmp(pRVE.Constantdirection,'Directions 2 and 3')
                REP2ndsize = sz(1);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = ( sz(1)/divisions(k) * sz(2) * sz(3) )^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).length = sz(1)/divisions(k);
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [sz(1)/divisions(k) sz(2) sz(3)] ./ sz(3);
                tmp_id = [];
                for x_=1:1:divisions(k)
                    x0 = floor((sz(1)/divisions(k)*(x_-1))+1);
                    x1 = floor(sz(1)/divisions(k)*x_);
                    y0 = 1;
                    y1 = sz(2);
                    z0 = 1;
                    z1 = sz(3);

                    sub_domain_id = sub_domain_id+1;
                    %actual_size = ((x1-x0+1)*(y1-y0+1)*(z1))^(1/3);
                    %actual_section = ((x1-x0+1)*(y1-y0+1))^(1/2);

                    All_subdomain(sub_domain_id,1) = sub_domain_id;
                    All_subdomain(sub_domain_id,2) = sub_domain_group;
                    All_subdomain(sub_domain_id,3) = number_subdomain_group;
                    All_subdomain(sub_domain_id,4) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length;
                    All_subdomain(sub_domain_id,5) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length;
                    All_subdomain(sub_domain_id,6) = GROUP_SUBDOMAIN.id(sub_domain_group).length;
                    All_subdomain(sub_domain_id,7) = x0;
                    All_subdomain(sub_domain_id,8) = x1;
                    All_subdomain(sub_domain_id,9) = y0;
                    All_subdomain(sub_domain_id,10) = y1;
                    All_subdomain(sub_domain_id,11) = z0;
                    All_subdomain(sub_domain_id,12) = z1;
                    All_subdomain(sub_domain_id,13) = divisions(k);

                    tmp_id = [tmp_id sub_domain_id];
                end
            end
            GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
        end
    end
end

if strcmp(pRVE.type,'E') || strcmp(pRVE.type,'F') % 'One subvolume + growing from volume center'

    if strcmp(pRVE.type,'E') % Keep initial aspect ratio
        d = 1/(pRVE.firstuniquevolume_size/100)^(1/dimension);
        center_volume_size = round(sz/d);
    elseif strcmp(pRVE.type,'F') % Custom aspect ratio
        v0 = prod(sz);
        v1 = pRVE.firstuniquevolume_size/100 * v0;
        k = (v1/prod(pRVE.Aspectratio_Onesub))^(1/dimension);
        center_volume_size = k*pRVE.Aspectratio_Onesub;
    end
    x0 = round(sz(1)/2 - center_volume_size(1)/2); x1 = round(sz(1)/2 + center_volume_size(1)/2);
    y0 = round(sz(2)/2 - center_volume_size(2)/2); y1 = round(sz(2)/2 + center_volume_size(2)/2);
    if dimension==3
        z0 = round(sz(3)/2 - center_volume_size(3)/2); z1 = round(sz(3)/2 + center_volume_size(3)/2);
    else
        z0 = 1; z1 = 1;
    end

    while 1
        if sub_domain_id>0
            g = pRVE.Growthperstep/100;
            if strcmp(pRVE.Growthrelativeto,'% of the previous subvolume')
                v0 = (x1-x0+1)*(y1-y0+1)*(z1-z0+1); % Current volume
                dv = v0*g;
            elseif strcmp(pRVE.Growthrelativeto,'% of the initial volume')
                dv = prod(sz)*g;
            end
            v0 = (x1-x0+1)*(y1-y0+1)*(z1-z0+1); % Current volume
            v1 = v0+dv; % New volume
            rv = v1/v0; % volume ratio
            rl = rv^(1/dimension); % lenght ratio
            lx = (x1-x0+1)*rl; ly = (y1-y0+1)*rl; lz = (z1-z0+1)*rl; % New length, preserve current aspect ratio
            x0 = round(sz(1)/2 - lx/2); x1 = round(sz(1)/2 + lx/2);
            y0 = round(sz(2)/2 - ly/2); y1 = round(sz(2)/2 + ly/2);
            if dimension == 3
                z0 = round(sz(3)/2 - lz/2); z1 = round(sz(3)/2 + lz/2);          
            else
                z0 = 1; z1 = 1;
            end
        end
        if x0<1 || y0<1 || z0<1 || x1>sz(1) || y1>sz(2) || z1>sz(3)
            break % Terminate while loop
        end
        
        sub_domain_group=sub_domain_group+1; % Increment group id
        number_subdomain_group=1; % Number of subdomain
        GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
        if dimension == 3
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = ( (x1-x0+1)*(y1-y0+1)*(z1-z0+1) ) ^(1/3);
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0; 
            GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [(x1-x0+1) (y1-y0+1) (z1-z0+1)]./(z1-z0+1);
        else
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = 0;
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = ( (x1-x0+1)*(y1-y0+1) ) ^(1/2); 
            GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [(x1-x0+1) (y1-y0+1)]./(y1-y0+1);
        end
        tmp_id = [];
        sub_domain_id = sub_domain_id+1;
        All_subdomain(sub_domain_id,1) = sub_domain_id;
        All_subdomain(sub_domain_id,2) = sub_domain_group;
        All_subdomain(sub_domain_id,3) = number_subdomain_group;
        All_subdomain(sub_domain_id,4) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length;
        All_subdomain(sub_domain_id,5) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length;
        All_subdomain(sub_domain_id,6) = GROUP_SUBDOMAIN.id(sub_domain_group).length;
        All_subdomain(sub_domain_id,7) = x0;
        All_subdomain(sub_domain_id,8) = x1;
        All_subdomain(sub_domain_id,9) = y0;
        All_subdomain(sub_domain_id,10) = y1;
        All_subdomain(sub_domain_id,11) = z0;
        All_subdomain(sub_domain_id,12) = z1;
        All_subdomain(sub_domain_id,13) = 0;
        tmp_id = [tmp_id sub_domain_id];
        GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
    end    
end

if strcmp(pRVE.type,'G') % 'One subvolume + growing from volume edge (G)'
    dinitial = round(sz * pRVE.firstuniquevolume_size/100);
    Representativity_2ndlength_str = 'length';
    if strcmp(pRVE.Growthdirection, 'Direction 1, from x min to x max')
        z=[1 sz(3)]; kz=[0 0]; kzz=0;
        y=[1 sz(2)]; ky=[0 0]; kyy=0;
        x=[1 dinitial(1)]; kx=[0 1]; kxx=1; 
        REP2ndsize = sz(1);        
    elseif strcmp(pRVE.Growthdirection,'Direction 1, from x max to x min')
        z=[1 sz(3)]; kz=[0 0]; kzz=0;
        y=[1 sz(2)]; ky=[0 0]; kyy=0;
        x=[sz(1)-dinitial(1)+1 sz(1)]; kx=[-1 0]; kxx=1; 
        REP2ndsize = sz(1);        
    elseif strcmp(pRVE.Growthdirection,'Direction 2, from y min to y max')
        z=[1 sz(3)]; kz=[0 0]; kzz=0;
        x=[1 sz(1)]; kx=[0 0]; kxx=0;
        y=[1 dinitial(2)]; ky=[0 1]; kyy=1; 
        REP2ndsize = sz(2);
    elseif strcmp(pRVE.Growthdirection,'Direction 2, from y max to y min')
        z=[1 sz(3)]; kz=[0 0]; kzz=0;
        x=[1 sz(1)]; kx=[0 0]; kxx=0;
        y=[sz(2)-dinitial(2)+1 sz(2)]; ky=[-1 0]; kyy=1;
        REP2ndsize = sz(2);        
    elseif strcmp(pRVE.Growthdirection,'Direction 3, from z min to z max')
        x=[1 sz(1)]; kx=[0 0]; kxx=0;
        y=[1 sz(2)]; ky=[0 0]; kyy=0;
        z=[1 dinitial(3)]; kz=[0 1]; kzz=1; 
        REP2ndsize = sz(3);        
    elseif strcmp(pRVE.Growthdirection,'Direction 3, from z max to z min')
        x=[1 sz(1)]; kx=[0 0]; kxx=0;
        y=[1 sz(2)]; ky=[0 0]; kyy=0;
        z=[sz(3)-dinitial(3)+1 sz(3)]; kz=[-1 0]; kzz=1; 
        REP2ndsize = sz(3);
    end
    x0=min(x); x1=max(x);
    y0=min(y); y1=max(y);
    if dimension == 3
        z0=min(z); z1=max(z);    
    else
        z0 = 1; z1 = 1;
    end
    
    while 1
        if sub_domain_id>0
            if strcmp(pRVE.Growthrelativeto,'% of the previous subvolume')
                current_subdomain_size = [x1-x0+1 y1-y0+1 z1-z0+1];
                dstep = round(pRVE.Growthperstep/100 * current_subdomain_size);
            elseif strcmp(pRVE.Growthrelativeto,'% of the initial volume')
                dstep = round(pRVE.Growthperstep/100 * sz);
            end
            x = x+(dstep(1)*kx);
            y = y+(dstep(2)*ky);
            z = z+(dstep(3)*kz);
        end
        x0=min(x); x1=max(x);
        y0=min(y); y1=max(y);
        if dimension == 3
            z0=min(z); z1=max(z);
        else
            z0 = 1; z1 = 1;
        end

        if x0<1 || y0<1 || z0<1 || x1>sz(1) || y1>sz(2) || z1>sz(3)
            break % Terminate while loop
        end
        
        sub_domain_group=sub_domain_group+1; % Increment group id
        number_subdomain_group=1; % Number of subdomain
        if dimension == 3
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = ( (x1-x0+1)*(y1-y0+1)*(z1-z0+1) ) ^(1/3);
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
            GROUP_SUBDOMAIN.id(sub_domain_group).length = abs( (x1-x0+1)*kxx + (y1-y0+1)*kyy + (z1-z0+1)*kzz); % Length;
            GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [(x1-x0+1) (y1-y0+1) (z1-z0+1)]./(z1-z0+1);
        else
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = 0;
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = ( (x1-x0+1)*(y1-y0+1) ) ^(1/2);
            GROUP_SUBDOMAIN.id(sub_domain_group).length = abs( (x1-x0+1)*kxx + (y1-y0+1)*kyy ); % Length;
            GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [(x1-x0+1) (y1-y0+1)]./(y1-y0+1);
        end
        tmp_id = [];
        sub_domain_id = sub_domain_id+1;
        All_subdomain(sub_domain_id,1) = sub_domain_id;
        All_subdomain(sub_domain_id,2) = sub_domain_group;
        All_subdomain(sub_domain_id,3) = number_subdomain_group;
        All_subdomain(sub_domain_id,4) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length;
        All_subdomain(sub_domain_id,5) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length;
        All_subdomain(sub_domain_id,6) = GROUP_SUBDOMAIN.id(sub_domain_group).length;
        All_subdomain(sub_domain_id,7) = x0;
        All_subdomain(sub_domain_id,8) = x1;
        All_subdomain(sub_domain_id,9) = y0;
        All_subdomain(sub_domain_id,10) = y1;
        All_subdomain(sub_domain_id,11) = z0;
        All_subdomain(sub_domain_id,12) = z1;
        All_subdomain(sub_domain_id,13) = 0;
        tmp_id = [tmp_id sub_domain_id];
        GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
    end
end


if strcmp(pRVE.type,'H') % 'One subvolume + growing with on constant length (H)'
    Constantdirection = pRVE.Growthdirection;
    d = (pRVE.firstuniquevolume_size/100)^(1/(dimension-1));
    if dimension == 3
        Representativity_2ndlength_str = 'square root';
    else
        Representativity_2ndlength_str = 'length';
    end

    if strcmp(Constantdirection,'Direction 1')
        REP2ndsize = (sz(2) * sz(3))^(1/(dimension-1));
        x0 = 1; x1 = sz(1);
        y0 = round(sz(2)/2 - d/2*sz(2)); y1 = round(sz(2)/2 + d/2*sz(2));
        z0 = round(sz(3)/2 - d/2*sz(3)); z1 = round(sz(3)/2 + d/2*sz(3));
    elseif strcmp(Constantdirection,'Direction 2')
        REP2ndsize = (sz(1) * sz(3))^(1/(dimension-1));
        x0 = round(sz(1)/2 - d/2*sz(1)); x1 = round(sz(1)/2 + d/2*sz(1));
        y0 = 1; y1 = sz(2);
        z0 = round(sz(3)/2 - d/2*sz(3)); z1 = round(sz(3)/2 + d/2*sz(3));
    elseif strcmp(Constantdirection,'Direction 3')
        REP2ndsize = (sz(1) * sz(2))^(1/(dimension-1));
        x0 = round(sz(1)/2 - d/2*sz(1)); x1 = round(sz(1)/2 + d/2*sz(1));
        y0 = round(sz(2)/2 - d/2*sz(2)); y1 = round(sz(2)/2 + d/2*sz(2));
        z0 = 1; z1 = sz(3);
    end

    if dimension == 2
        z0 = 1; z1 = 1;
    end    
    
    while 1
        if sub_domain_id>0
            g = pRVE.Growthperstep/100;
            if strcmp(pRVE.Growthrelativeto,'% of the previous subvolume')
                v0 = (x1-x0+1)*(y1-y0+1)*(z1-z0+1); % Current volume
                dv = v0*g;
            elseif strcmp(pRVE.Growthrelativeto,'% of the initial volume')
                dv = prod(sz)*g;
            end 
            v0 = (x1-x0+1)*(y1-y0+1)*(z1-z0+1); % Current volume
            v1 = v0+dv; % New volume
            d = (v1/prod(sz))^(1/(dimension-1));
            if strcmp(Constantdirection,'Direction 1')
                x0 = 1; x1 = sz(1);
                y0 = round(sz(2)/2 - d/2*sz(2)); y1 = round(sz(2)/2 + d/2*sz(2));
                z0 = round(sz(3)/2 - d/2*sz(3)); z1 = round(sz(3)/2 + d/2*sz(3));
            elseif strcmp(Constantdirection,'Direction 2')
                x0 = round(sz(1)/2 - d/2*sz(1)); x1 = round(sz(1)/2 + d/2*sz(1));
                y0 = 1; y1 = sz(2);
                z0 = round(sz(3)/2 - d/2*sz(3)); z1 = round(sz(3)/2 + d/2*sz(3));
            elseif strcmp(Constantdirection,'Direction 3')
                x0 = round(sz(1)/2 - d/2*sz(1)); x1 = round(sz(1)/2 + d/2*sz(1));
                y0 = round(sz(2)/2 - d/2*sz(2)); y1 = round(sz(2)/2 + d/2*sz(2));
                z0 = 1; z1 = sz(3);
            end
            if dimension == 2
                z0 = 1; z1 = 1;
            end
        end
        if x0<1 || y0<1 || z0<1 || x1>sz(1) || y1>sz(2) || z1>sz(3)
            break % Terminate while loop
        end

        sub_domain_group=sub_domain_group+1; % Increment group id
        number_subdomain_group=1; % Number of subdomain

        if dimension == 3
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = ( (x1-x0+1)*(y1-y0+1)*(z1-z0+1) ) ^(1/3);
            if strcmp(Constantdirection,'Direction 1')
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = ((y1-y0+1)*(z1-z0+1))^(1/2);
            elseif strcmp(Constantdirection,'Direction 2')
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = ((x1-x0+1)*(z1-z0+1))^(1/2);
            elseif strcmp(Constantdirection,'Direction 3')
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = ((x1-x0+1)*(y1-y0+1))^(1/2);
            end
            GROUP_SUBDOMAIN.id(sub_domain_group).length = 0; % Length;
            GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [(x1-x0+1) (y1-y0+1) (z1-z0+1)]./(z1-z0+1);
        else
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = 0;
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = ((x1-x0+1)*(y1-y0+1))^(1/2);
            if strcmp(Constantdirection,'Direction 1')
                GROUP_SUBDOMAIN.id(sub_domain_group).length = (y1-y0+1);
            elseif strcmp(Constantdirection,'Direction 2')
                GROUP_SUBDOMAIN.id(sub_domain_group).length = (x1-x0+1);
            elseif strcmp(Constantdirection,'Direction 3')
                GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
            end
            GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [(x1-x0+1) (y1-y0+1)]./(y1-y0+1);
        end

        tmp_id = [];
        sub_domain_id = sub_domain_id+1;
        All_subdomain(sub_domain_id,1) = sub_domain_id;
        All_subdomain(sub_domain_id,2) = sub_domain_group;
        All_subdomain(sub_domain_id,3) = number_subdomain_group;
        All_subdomain(sub_domain_id,4) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length;
        All_subdomain(sub_domain_id,5) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length;
        All_subdomain(sub_domain_id,6) = GROUP_SUBDOMAIN.id(sub_domain_group).length;
        All_subdomain(sub_domain_id,7) = x0;
        All_subdomain(sub_domain_id,8) = x1;
        All_subdomain(sub_domain_id,9) = y0;
        All_subdomain(sub_domain_id,10) = y1;
        All_subdomain(sub_domain_id,11) = z0;
        All_subdomain(sub_domain_id,12) = z1;
        All_subdomain(sub_domain_id,13) = 0;
        tmp_id = [tmp_id sub_domain_id];
        GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
    end
end

if ~isempty(Representativity_2ndlength_str)
    REP2ndsize = REP2ndsize * voxel_size;
    if strcmp(Representativity_2ndlength_str,'square root')
        if strcmp(pRVE.analysis,'Independent subvolumes')
            Representativity_analysis_size2 = 'Representative Section Area (RSA)';
        else
            Representativity_analysis_size2 = 'Section area convergence';
        end
    elseif strcmp(Representativity_2ndlength_str,'length')
        if strcmp(pRVE.analysis,'Independent subvolumes')
            Representativity_analysis_size2 = 'Representative length (RL)';
        else
            Representativity_analysis_size2 = 'Length convergence';
        end
    end
end

end




