function [All_subdomain,GROUP_SUBDOMAIN, Wholevolume] = Function_get_Subdomains(RVEparameters, Domain_size, Initial_Domain_size)
% Determine the bounds of subdomains from a domain.

sub_domain_id=0;
sub_domain_group=0;

divisions = RVEparameters.divisions;

Wholevolume(1) = (prod(Domain_size))^(1/3);
Wholevolume(3) = Domain_size(1)/Domain_size(3); % Aspect ratio
Wholevolume(4) = Domain_size(2)/Domain_size(3);
Wholevolume(5) = Domain_size(3)/Domain_size(3);

if strcmp(RVEparameters.type,'A') % 'Independant subvolumes of same size + keep initial aspect ratio (A)'

    initial_volume = prod(Initial_Domain_size);
    current_volume = prod(Domain_size);
    minimum_subvolumesize = RVEparameters.minsubsize/100 * initial_volume;
    max_division = round((current_volume/minimum_subvolumesize)^(1/3));
    if strcmp(RVEparameters.iterateuntil,'Minimum subvolume size has been reached (ignore divider value)')
        divisions = 2:1:max_division;
    else
        divisions(divisions>max_division) = [];
    end

    if ~isempty(divisions)
        for k=1:1:length(divisions)
            sub_domain_group=sub_domain_group+1;
            number_subdomain_group=divisions(k)^3;
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1)/divisions(k) * Domain_size(2)/divisions(k) * Domain_size(3)/divisions(k))^(1/3);
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
            GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
            GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1)/divisions(k) Domain_size(2)/divisions(k) Domain_size(3)/divisions(k)] ./ (Domain_size(3)/divisions(k));
            tmp_id = [];
            for x_=1:1:divisions(k)
                for y_=1:1:divisions(k)
                    for z_=1:1:divisions(k)

                        x0 = floor((Domain_size(1)/divisions(k)*(x_-1))+1);
                        x1 = floor(Domain_size(1)/divisions(k)*x_);
                        y0 = floor((Domain_size(2)/divisions(k)*(y_-1))+1);
                        y1 = floor(Domain_size(2)/divisions(k)*y_);
                        z0 = floor((Domain_size(3)/divisions(k)*(z_-1))+1);
                        z1 = floor(Domain_size(3)/divisions(k)*z_);

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
                        All_subdomain(sub_domain_id,13) = divisions(k);

                        tmp_id = [tmp_id sub_domain_id];
                    end
                end
            end
            GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
        end
    end
    
    if RVEparameters.subs4
        if current_volume/4 > minimum_subvolumesize
            sub_domain_group=sub_domain_group+1; % Increment group id
            number_subdomain_group=4; % Number of subdomain
            tmp_id = [];
            if Domain_size(1)<=Domain_size(2) && Domain_size(1)<=Domain_size(3)
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(2)/2 * Domain_size(3)/2 * Domain_size(1))^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1) Domain_size(2)/2 Domain_size(3)/2] ./ (Domain_size(3)/2);
                nx = 1; ny=2; nz=2;
                yy_start(1)=1; yy_end(1)=floor(Domain_size(2)/2);
                yy_start(2)=yy_end(1)+1; yy_end(2)=Domain_size(2);
                zz_start(1)=1; zz_end(1)=floor(Domain_size(3)/2);
                zz_start(2)=zz_end(1)+1; zz_end(2)=Domain_size(3);
                xx_start(1) = 1; xx_end(1) = Domain_size(1);
            end

            if Domain_size(2)<=Domain_size(1) && Domain_size(2)<=Domain_size(3)
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1)/2 * Domain_size(3)/2 * Domain_size(2))^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1)/2 Domain_size(2) Domain_size(3)/2] ./ (Domain_size(3)/2);
                nx = 2; ny=1; nz=2;
                xx_start(1)=1; xx_end(1)=floor(Domain_size(1)/2);
                xx_start(2)=xx_end(1)+1; xx_end(2)=Domain_size(1);
                zz_start(1)=1; zz_end(1)=floor(Domain_size(3)/2);
                zz_start(2)=zz_end(1)+1; zz_end(2)=Domain_size(3);
                yy_start(1) = 1; yy_end(1) = Domain_size(2);
            end

            if Domain_size(3)<=Domain_size(1) && Domain_size(3)<=Domain_size(2)
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1)/2 * Domain_size(2)/2 * Domain_size(3))^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1)/2 Domain_size(2)/2 Domain_size(3)] ./ (Domain_size(3));
                nx = 2; ny=2; nz=1;
                xx_start(1)=1; xx_end(1)=floor(Domain_size(1)/2);
                xx_start(2)=xx_end(1)+1; xx_end(2)=Domain_size(1);
                yy_start(1)=1; yy_end(1)=floor(Domain_size(2)/2);
                yy_start(2)=yy_end(1)+1; yy_end(2)=Domain_size(2);
                zz_start(1) = 1; zz_end(1) = Domain_size(3);
            end

            for kx=1:1:nx
                for ky=1:1:ny
                    for kz=1:1:nz
                        sub_domain_id = sub_domain_id+1;
                        All_subdomain(sub_domain_id,1) = sub_domain_id;
                        All_subdomain(sub_domain_id,2) = sub_domain_group;
                        All_subdomain(sub_domain_id,3) = number_subdomain_group;
                        All_subdomain(sub_domain_id,4) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length;
                        All_subdomain(sub_domain_id,5) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length;
                        All_subdomain(sub_domain_id,6) = GROUP_SUBDOMAIN.id(sub_domain_group).length;
                        All_subdomain(sub_domain_id,7) = xx_start(kx);
                        All_subdomain(sub_domain_id,8) = xx_end(kx);
                        All_subdomain(sub_domain_id,9) = yy_start(ky);
                        All_subdomain(sub_domain_id,10) = yy_end(ky);
                        All_subdomain(sub_domain_id,11) = zz_start(kz);
                        All_subdomain(sub_domain_id,12) = zz_end(kz);
                        All_subdomain(sub_domain_id,13) = 0;
                        tmp_id = [tmp_id sub_domain_id];
                    end
                end
            end

            GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
        end
    end

    if RVEparameters.subs2
        if current_volume/2 > minimum_subvolumesize
            sub_domain_group=sub_domain_group+1; % Increment group id
            number_subdomain_group=2; % Number of subdomain
            tmp_id = [];
            if Domain_size(1)<=Domain_size(2) && Domain_size(1)<=Domain_size(3)
                if Domain_size(2)>=Domain_size(3)
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(2)/2 * Domain_size(3) * Domain_size(1))^(1/3);
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1) Domain_size(2)/2 Domain_size(3)] ./ Domain_size(3);
                    yy_start(1)=1; yy_end(1)=floor(Domain_size(2)/2);
                    yy_start(2)=yy_end(1)+1; yy_end(2)=Domain_size(2);
                    zz_start(1)=1; zz_end(1)=Domain_size(3);
                    zz_start(2)=1; zz_end(2)=Domain_size(3);
                else
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(2) * Domain_size(3)/2 * Domain_size(1))^(1/3);
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1) Domain_size(2) Domain_size(3)/2] ./ (Domain_size(3)/2);
                    yy_start(1)=1; yy_end(1)=Domain_size(2);
                    yy_start(2)=1; yy_end(2)=Domain_size(2);
                    zz_start(1)=1; zz_end(1)=floor(Domain_size(3)/2);
                    zz_start(2)=zz_end(1)+1; zz_end(2)=Domain_size(3);
                end
                xx_start(1) = 1; xx_end(1) = Domain_size(1);
                xx_start(2) = 1; xx_end(2) = Domain_size(1);
            end

            if Domain_size(2)<=Domain_size(1) && Domain_size(2)<=Domain_size(3)
                if Domain_size(1)>=Domain_size(3)
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1)/2 * Domain_size(3) * Domain_size(2))^(1/3);
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1)/2 Domain_size(2) Domain_size(3)] ./ Domain_size(3);
                    xx_start(1)=1; xx_end(1)=floor(Domain_size(1)/2);
                    xx_start(2)=xx_end(1)+1; xx_end(2)=Domain_size(1);
                    zz_start(1)=1; zz_end(1)=Domain_size(3);
                    zz_start(2)=1; zz_end(2)=Domain_size(3);
                else
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1) * Domain_size(3)/2 * Domain_size(2))^(1/3);
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1) Domain_size(2) Domain_size(3)/2] ./ (Domain_size(3)/2);
                    xx_start(1)=1; xx_end(1)=Domain_size(1);
                    xx_start(2)=1; xx_end(2)=Domain_size(1);
                    zz_start(1)=1; zz_end(1)=floor(Domain_size(3)/2);
                    zz_start(2)=zz_end(1)+1; zz_end(2)=Domain_size(3);
                end
                yy_start(1) = 1; yy_end(1) = Domain_size(2);
                yy_start(2) = 1; yy_end(2) = Domain_size(2);
            end

            if Domain_size(3)<=Domain_size(1) && Domain_size(3)<=Domain_size(2)
                if Domain_size(1)>=Domain_size(2)
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1)/2 * Domain_size(2) * Domain_size(3))^(1/3);
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1)/2 Domain_size(2) Domain_size(3)] ./ Domain_size(3);
                    xx_start(1)=1; xx_end(1)=floor(Domain_size(1)/2);
                    xx_start(2)=xx_end(1)+1; xx_end(2)=Domain_size(1);
                    yy_start(1)=1; yy_end(1)=Domain_size(2);
                    yy_start(2)=1; yy_end(2)=Domain_size(2);
                else
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1) * Domain_size(2)/2 * Domain_size(3))^(1/3);
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1) Domain_size(2)/2 Domain_size(3)] ./ Domain_size(3);
                    xx_start(1)=1; xx_end(1)=Domain_size(1);
                    xx_start(2)=1; xx_end(2)=Domain_size(1);
                    yy_start(1)=1; yy_end(1)=floor(Domain_size(2)/2);
                    yy_start(2)=yy_end(1)+1; yy_end(2)=Domain_size(2);
                end
                zz_start(1) = 1; zz_end(1) = Domain_size(3);
                zz_start(2) = 1; zz_end(2) = Domain_size(3);
            end

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


if  strcmp(RVEparameters.type,'B') % 'Independant subvolumes of same size + user-defined aspect ratio (B)'
    
    AR_desired = RVEparameters.Aspectratio / RVEparameters.Aspectratio(3);
    
    % Search for the largest subvolume that has the desired AR
    l1_n = 1; l1(1) = Initial_Domain_size(1)/l1_n;
    l2(1) = AR_desired(2)*l1(1)/AR_desired(1); l2_n = Initial_Domain_size(2)/l2(1);
    l3(1) = AR_desired(3)*l1(1)/AR_desired(1); l3_n = Initial_Domain_size(3)/l3(1);
    vol = [0 0 0];
    if l2_n>=1 && l3_n>=1
        vol(1)=prod([l1(1) l2(1) l3(1)]);
    end

    l2_n = 1; l2(2) = Initial_Domain_size(2)/l2_n;
    l1(2) = AR_desired(1)*l2(2)/AR_desired(2); l1_n = Initial_Domain_size(1)/l1(2);
    l3(2) = AR_desired(3)*l2(2)/AR_desired(2); l3_n = Initial_Domain_size(3)/l3(2);
    if l1_n>=1 && l3_n>=1
        vol(2)=prod([l1(2) l2(2) l3(2)]);
    end

    l3_n = 1; l3(3) = Initial_Domain_size(3)/l3_n;
    l1(3) = AR_desired(1)*l3(3)/AR_desired(3); l1_n = Initial_Domain_size(1)/l1(3);
    l2(3) = AR_desired(2)*l3(3)/AR_desired(3); l2_n = Initial_Domain_size(2)/l2(3);
    if l1_n>=1 && l2_n>=1
        vol(3)=prod([l1(3) l2(3) l3(3)]);
    end

    idx=find(vol==max(vol)); % Case that ensure the largest volume
    idx=idx(1); % In case of equality
    Max_subdomain_size = [l1(idx) l2(idx) l3(idx)];

    initial_volume = prod(Max_subdomain_size);
    minimum_subvolumesize = RVEparameters.minsubsize/100 * initial_volume;
    max_division = round((initial_volume/minimum_subvolumesize)^(1/3));

    if sum(Domain_size~=Initial_Domain_size)
        % Search for the largest subvolume that has the desired AR
        l1_n = 1; l1(1) = Domain_size(1)/l1_n;
        l2(1) = AR_desired(2)*l1(1)/AR_desired(1); l2_n = Domain_size(2)/l2(1);
        l3(1) = AR_desired(3)*l1(1)/AR_desired(1); l3_n = Domain_size(3)/l3(1);
        vol = [0 0 0];
        if l2_n>=1 && l3_n>=1
            vol(1)=prod([l1(1) l2(1) l3(1)]);
        end

        l2_n = 1; l2(2) = Domain_size(2)/l2_n;
        l1(2) = AR_desired(1)*l2(2)/AR_desired(2); l1_n = Domain_size(1)/l1(2);
        l3(2) = AR_desired(3)*l2(2)/AR_desired(2); l3_n = Domain_size(3)/l3(2);
        if l1_n>=1 && l3_n>=1
            vol(2)=prod([l1(2) l2(2) l3(2)]);
        end

        l3_n = 1; l3(3) = Domain_size(3)/l3_n;
        l1(3) = AR_desired(1)*l3(3)/AR_desired(3); l1_n = Domain_size(1)/l1(3);
        l2(3) = AR_desired(2)*l3(3)/AR_desired(3); l2_n = Domain_size(2)/l2(3);
        if l1_n>=1 && l2_n>=1
            vol(3)=prod([l1(3) l2(3) l3(3)]);
        end

        idx=find(vol==max(vol)); % Case that ensure the largest volume
        idx=idx(1); % In case of equality
        Max_subdomain_size = [l1(idx) l2(idx) l3(idx)];
        max_division = round( (prod(Max_subdomain_size) /minimum_subvolumesize)^(1/3));
    end

    if strcmp(RVEparameters.iterateuntil,'Minimum subvolume size has been reached (ignore divider value)')
        divisions = 2:1:max_division;
    else
        divisions(divisions>max_division) = [];
    end

    for k=0:1:length(divisions)
        sub_domain_group=sub_domain_group+1;
        if k==0
            current_subdomain_size = Max_subdomain_size;
        else
            current_subdomain_size = Max_subdomain_size/divisions(k);
        end
        div_ = floor(Domain_size./current_subdomain_size);
        current_domain_size = div_.*current_subdomain_size;
        number_subdomain_group=prod(div_);
        
        GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (current_domain_size(1)/div_(1) * current_domain_size(2)/div_(2) * current_domain_size(3)/div_(3))^(1/3);
        GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
        GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
        GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [current_domain_size(1)/div_(1) current_domain_size(2)/div_(2) current_domain_size(3)/div_(3)] ./ (current_domain_size(3)/div_(3));
        tmp_id = [];
        for x_=1:1:div_(1)
            for y_=1:1:div_(2)
                for z_=1:1:div_(3)

                    x0 = floor((current_domain_size(1)/div_(1)*(x_-1))+1);
                    x1 = floor(current_domain_size(1)/div_(1)*x_);
                    y0 = floor((current_domain_size(2)/div_(2)*(y_-1))+1);
                    y1 = floor(current_domain_size(2)/div_(2)*y_);
                    z0 = floor((current_domain_size(3)/div_(3)*(z_-1))+1);
                    z1 = floor(current_domain_size(3)/div_(3)*z_);
                    
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


if strcmp(RVEparameters.type,'C') % 'Independant subvolumes of same size + one constant length (C)'
    
    initial_volume = prod(Initial_Domain_size);
    current_volume = prod(Domain_size);
    minimum_subvolumesize = RVEparameters.minsubsize/100 * initial_volume;
    max_division = round((current_volume/minimum_subvolumesize)^(1/2));
    if strcmp(RVEparameters.iterateuntil,'Minimum subvolume size has been reached (ignore divider value)')
        divisions = 2:1:max_division;
    else
        divisions(divisions>max_division) = [];
    end    

    if ~isempty(divisions)
        for k=1:1:length(divisions)
            sub_domain_group=sub_domain_group+1; % Increment group id
            number_subdomain_group=divisions(k)^2; % Number of subdomain
            if strcmp(RVEparameters.Constantdirection,'Direction 1')
                Wholevolume(2) = (Domain_size(2) * Domain_size(3))^(1/2);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(2)/divisions(k) * Domain_size(3)/divisions(k) * Domain_size(1))^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (Domain_size(2)/divisions(k) * Domain_size(3)/divisions(k))^(1/2);
                GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1) Domain_size(2)/divisions(k) Domain_size(3)/divisions(k)] ./ (Domain_size(3)/divisions(k));
                tmp_id = [];
                for y_=1:1:divisions(k)
                    for z_=1:1:divisions(k)
                        x0 = 1;
                        x1 = Domain_size(1);
                        y0 = floor((Domain_size(2)/divisions(k)*(y_-1))+1);
                        y1 = floor(Domain_size(2)/divisions(k)*y_);
                        z0 = floor((Domain_size(3)/divisions(k)*(z_-1))+1);
                        z1 = floor(Domain_size(3)/divisions(k)*z_);

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

            elseif strcmp(RVEparameters.Constantdirection,'Direction 2')
                Wholevolume(2) = (Domain_size(1) * Domain_size(3))^(1/2);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1)/divisions(k) * Domain_size(3)/divisions(k) * Domain_size(2))^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (Domain_size(1)/divisions(k) * Domain_size(3)/divisions(k))^(1/2);
                GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1)/divisions(k) Domain_size(2) Domain_size(3)/divisions(k)] ./ (Domain_size(3)/divisions(k));
                tmp_id = [];
                for x_=1:1:divisions(k)
                    for z_=1:1:divisions(k)
                        x0 = floor((Domain_size(1)/divisions(k)*(x_-1))+1);
                        x1 = floor(Domain_size(1)/divisions(k)*x_);
                        y0 = 1;
                        y1 = Domain_size(2);
                        z0 = floor((Domain_size(3)/divisions(k)*(z_-1))+1);
                        z1 = floor(Domain_size(3)/divisions(k)*z_);

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

            elseif strcmp(RVEparameters.Constantdirection,'Direction 3')
                Wholevolume(2) = (Domain_size(1) * Domain_size(2))^(1/2);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1)/divisions(k) * Domain_size(2)/divisions(k) * Domain_size(3))^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (Domain_size(1)/divisions(k) * Domain_size(2)/divisions(k))^(1/2);
                GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1)/divisions(k) Domain_size(2)/divisions(k) Domain_size(3)] ./ Domain_size(3);
                tmp_id = [];
                for x_=1:1:divisions(k)
                    for y_=1:1:divisions(k)
                        x0 = floor((Domain_size(1)/divisions(k)*(x_-1))+1);
                        x1 = floor(Domain_size(1)/divisions(k)*x_);
                        y0 = floor((Domain_size(2)/divisions(k)*(y_-1))+1);
                        y1 = floor(Domain_size(2)/divisions(k)*y_);
                        z0 = 1;
                        z1 = Domain_size(3);

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
            end
            GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
        end
    end

    if RVEparameters.subs2
        if current_volume/2 > minimum_subvolumesize
            sub_domain_group=sub_domain_group+1; % Increment group id
            number_subdomain_group=2; % Number of subdomain
            tmp_id = [];
            if strcmp(RVEparameters.Constantdirection,'Direction 1')
                if Domain_size(2)>=Domain_size(3)
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(2)/2 * Domain_size(3) * Domain_size(1))^(1/3);
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (Domain_size(2)/2 * Domain_size(3))^(1/2);
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1) Domain_size(2)/2 Domain_size(3)] ./ Domain_size(3);
                    yy_start(1)=1; yy_end(1)=floor(Domain_size(2)/2);
                    yy_start(2)=yy_end(1)+1; yy_end(2)=Domain_size(2);
                    zz_start(1)=1; zz_end(1)=Domain_size(3);
                    zz_start(2)=1; zz_end(2)=Domain_size(3);
                else
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(2) * Domain_size(3)/2 * Domain_size(1))^(1/3);
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (Domain_size(2) * Domain_size(3)/2)^(1/2);
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1) Domain_size(2) Domain_size(3)/2] ./ (Domain_size(3)/2);
                    yy_start(1)=1; yy_end(1)=Domain_size(2);
                    yy_start(2)=1; yy_end(2)=Domain_size(2);
                    zz_start(1)=1; zz_end(1)=floor(Domain_size(3)/2);
                    zz_start(2)=zz_end(1)+1; zz_end(2)=Domain_size(3);
                end
                xx_start(1) = 1; xx_end(1) = Domain_size(1);
                xx_start(2) = 1; xx_end(2) = Domain_size(1);

            elseif strcmp(RVEparameters.Constantdirection,'Direction 2')
                if Domain_size(1)>=Domain_size(3)
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1)/2 * Domain_size(3) * Domain_size(2))^(1/3);
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (Domain_size(1)/2 * Domain_size(3))^(1/2);
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1)/2 Domain_size(2) Domain_size(3)] ./ Domain_size(3);
                    xx_start(1)=1; xx_end(1)=floor(Domain_size(1)/2);
                    xx_start(2)=xx_end(1)+1; xx_end(2)=Domain_size(1);
                    zz_start(1)=1; zz_end(1)=Domain_size(3);
                    zz_start(2)=1; zz_end(2)=Domain_size(3);
                else
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1) * Domain_size(3)/2 * Domain_size(2))^(1/3);
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (Domain_size(1) * Domain_size(3)/2)^(1/2);
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1) Domain_size(2) Domain_size(3)/2] ./ (Domain_size(3)/2);
                    xx_start(1)=1; xx_end(1)=Domain_size(1);
                    xx_start(2)=1; xx_end(2)=Domain_size(1);
                    zz_start(1)=1; zz_end(1)=floor(Domain_size(3)/2);
                    zz_start(2)=zz_end(1)+1; zz_end(2)=Domain_size(3);
                end
                yy_start(1) = 1; yy_end(1) = Domain_size(2);
                yy_start(2) = 1; yy_end(2) = Domain_size(2);

            elseif strcmp(RVEparameters.Constantdirection,'Direction 3')
                if Domain_size(1)>=Domain_size(2)
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1)/2 * Domain_size(2) * Domain_size(3))^(1/3);
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (Domain_size(1)/2 * Domain_size(2))^(1/2);
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1)/2 Domain_size(2) Domain_size(3)] ./ Domain_size(3);
                    xx_start(1)=1; xx_end(1)=floor(Domain_size(1)/2);
                    xx_start(2)=xx_end(1)+1; xx_end(2)=Domain_size(1);
                    yy_start(1)=1; yy_end(1)=Domain_size(2);
                    yy_start(2)=1; yy_end(2)=Domain_size(2);
                else
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1) * Domain_size(2)/2 * Domain_size(3))^(1/3);
                    GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (Domain_size(1) * Domain_size(2)/2)^(1/2);
                    GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
                    GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1) Domain_size(2)/2 Domain_size(3)] ./ Domain_size(3);
                    xx_start(1)=1; xx_end(1)=Domain_size(1);
                    xx_start(2)=1; xx_end(2)=Domain_size(1);
                    yy_start(1)=1; yy_end(1)=floor(Domain_size(2)/2);
                    yy_start(2)=yy_end(1)+1; yy_end(2)=Domain_size(2);
                end
                zz_start(1) = 1; zz_end(1) = Domain_size(3);
                zz_start(2) = 1; zz_end(2) = Domain_size(3);
            end
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

if strcmp(RVEparameters.type,'D') % 'Independant subvolumes of same size + two constant lengths (D)'
    
    initial_volume = prod(Initial_Domain_size);
    current_volume = prod(Domain_size);
    minimum_subvolumesize = RVEparameters.minsubsize/100 * initial_volume;
    max_division = round((current_volume/minimum_subvolumesize)^(1/1));
    if strcmp(RVEparameters.iterateuntil,'Minimum subvolume size has been reached (ignore divider value)')
        divisions = 2:1:max_division;
    else
        divisions(divisions>max_division) = [];
    end    

    if ~isempty(divisions)
        for k=1:1:length(divisions)
            sub_domain_group=sub_domain_group+1; % Increment group id
            number_subdomain_group=divisions(k)^1; % Number of subdomain
            if strcmp(RVEparameters.Constantdirection,'Direction 1 and 2')
                Wholevolume(2) = Domain_size(3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = ( Domain_size(1) * Domain_size(2) * Domain_size(3)/divisions(k) )^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).length = Domain_size(3)/divisions(k);
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1) Domain_size(2) Domain_size(3)/divisions(k)] ./ (Domain_size(3)/divisions(k));
                tmp_id = [];
                for z_=1:1:divisions(k)
                    x0 = 1;
                    x1 = Domain_size(1);
                    y0 = 1;
                    y1 = Domain_size(2);
                    z0 = floor((Domain_size(3)/divisions(k)*(z_-1))+1);
                    z1 = floor(Domain_size(3)/divisions(k)*z_);

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

            elseif strcmp(RVEparameters.Constantdirection,'Direction 1 and 3')
                Wholevolume(2) = Domain_size(2);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = ( Domain_size(1) * Domain_size(2)/divisions(k) * Domain_size(3) )^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).length = Domain_size(2)/divisions(k);
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1) Domain_size(2)/divisions(k) Domain_size(3)] ./ Domain_size(3);
                tmp_id = [];
                for y_=1:1:divisions(k)
                    x0 = 1;
                    x1 = Domain_size(1);
                    y0 = floor((Domain_size(2)/divisions(k)*(y_-1))+1);
                    y1 = floor(Domain_size(2)/divisions(k)*y_);
                    z0 = 1;
                    z1 = Domain_size(3);

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

            elseif strcmp(RVEparameters.Constantdirection,'Direction 2 and 3')
                Wholevolume(2) = Domain_size(1);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = ( Domain_size(1)/divisions(k) * Domain_size(2) * Domain_size(3) )^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).length = Domain_size(1)/divisions(k);
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1)/divisions(k) Domain_size(2) Domain_size(3)] ./ Domain_size(3);
                tmp_id = [];
                for x_=1:1:divisions(k)
                    x0 = floor((Domain_size(1)/divisions(k)*(x_-1))+1);
                    x1 = floor(Domain_size(1)/divisions(k)*x_);
                    y0 = 1;
                    y1 = Domain_size(2);
                    z0 = 1;
                    z1 = Domain_size(3);

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

if strcmp(RVEparameters.type,'E') || strcmp(RVEparameters.type,'F') % 'One subvolume + growing from volume center'
    if strcmp(RVEparameters.type,'E') % Keep initial aspect ratio
        d = 1/(RVEparameters.firstuniquevolume_size/100)^(1/3);
        center_volume_size = round(Domain_size/d);
    elseif strcmp(RVEparameters.type,'F') % Custom aspect ratio
        v0 = prod(Domain_size);
        v1 = RVEparameters.firstuniquevolume_size/100 * v0;
        k = (v1/prod(RVEparameters.Aspectratio))^(1/3);
        center_volume_size = k*RVEparameters.Aspectratio;
    end
    x0 = round(Domain_size(1)/2 - center_volume_size(1)/2); x1 = round(Domain_size(1)/2 + center_volume_size(1)/2);
    y0 = round(Domain_size(2)/2 - center_volume_size(2)/2); y1 = round(Domain_size(2)/2 + center_volume_size(2)/2);
    z0 = round(Domain_size(3)/2 - center_volume_size(3)/2); z1 = round(Domain_size(3)/2 + center_volume_size(3)/2);

    while 1
        if sub_domain_id>0
            g = RVEparameters.Growthperstep/100;
            if strcmp(RVEparameters.Growthrelativeto,'% of the previous subvolume')
                v0 = (x1-x0+1)*(y1-y0+1)*(z1-z0+1); % Current volume
                dv = v0*g;
            elseif strcmp(RVEparameters.Growthrelativeto,'% of the full volume')
                dv = prod(Domain_size)*g;
            end
            v0 = (x1-x0+1)*(y1-y0+1)*(z1-z0+1); % Current volume
            v1 = v0+dv; % New volume
            rv = v1/v0; % volume ratio
            rl = rv^(1/3); % lenght ratio
            lx = (x1-x0+1)*rl; ly = (y1-y0+1)*rl; lz = (z1-z0+1)*rl; % New length, preserve current aspect ratio
            x0 = round(Domain_size(1)/2 - lx/2); x1 = round(Domain_size(1)/2 + lx/2);
            y0 = round(Domain_size(2)/2 - ly/2); y1 = round(Domain_size(2)/2 + ly/2);
            z0 = round(Domain_size(3)/2 - lz/2); z1 = round(Domain_size(3)/2 + lz/2);          
        end
        if x0<1 || y0<1 || z0<1 || x1>Domain_size(1) || y1>Domain_size(2) || z1>Domain_size(3)
            break % Terminate while loop
        end
        
        sub_domain_group=sub_domain_group+1; % Increment group id
        number_subdomain_group=1; % Number of subdomain
        GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = ( (x1-x0+1)*(y1-y0+1)*(z1-z0+1) ) ^(1/3);
        GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0; % Length
        GROUP_SUBDOMAIN.id(sub_domain_group).length = 0;
        GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [(x1-x0+1) (y1-y0+1) (z1-z0+1)]./(z1-z0+1);
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

if strcmp(RVEparameters.type,'G') % 'One subvolume + growing from volume edge (G)'
    dinitial = round(Domain_size * RVEparameters.firstuniquevolume_size/100);
    if strcmp(RVEparameters.Growthdirection, 'Direction 1, from x min to x max')
        z=[1 Domain_size(3)]; kz=[0 0]; kzz=0;
        y=[1 Domain_size(2)]; ky=[0 0]; kyy=0;
        x=[1 dinitial(1)]; kx=[0 1]; kxx=1; 
        Wholevolume(2) = Domain_size(1);        
    elseif strcmp(RVEparameters.Growthdirection,'Direction 1, from x max to x min')
        z=[1 Domain_size(3)]; kz=[0 0]; kzz=0;
        y=[1 Domain_size(2)]; ky=[0 0]; kyy=0;
        x=[Domain_size(1)-dinitial(1)+1 Domain_size(1)]; kx=[-1 0]; kxx=1; 
        Wholevolume(2) = Domain_size(1);        
    elseif strcmp(RVEparameters.Growthdirection,'Direction 2, from y min to y max')
        z=[1 Domain_size(3)]; kz=[0 0]; kzz=0;
        x=[1 Domain_size(1)]; kx=[0 0]; kxx=0;
        y=[1 dinitial(2)]; ky=[0 1]; kyy=1; 
        Wholevolume(2) = Domain_size(2);
    elseif strcmp(RVEparameters.Growthdirection,'Direction 2, from y max to y min')
        z=[1 Domain_size(3)]; kz=[0 0]; kzz=0;
        x=[1 Domain_size(1)]; kx=[0 0]; kxx=0;
        y=[Domain_size(2)-dinitial(2)+1 Domain_size(2)]; ky=[-1 0]; kyy=1;
        Wholevolume(2) = Domain_size(2);        
    elseif strcmp(RVEparameters.Growthdirection,'Direction 3, from z min to z max')
        x=[1 Domain_size(1)]; kx=[0 0]; kxx=0;
        y=[1 Domain_size(2)]; ky=[0 0]; kyy=0;
        z=[1 dinitial(3)]; kz=[0 1]; kzz=1; 
        Wholevolume(2) = Domain_size(3);        
    elseif strcmp(RVEparameters.Growthdirection,'Direction 3, from z max to z min')
        x=[1 Domain_size(1)]; kx=[0 0]; kxx=0;
        y=[1 Domain_size(2)]; ky=[0 0]; kyy=0;
        z=[Domain_size(3)-dinitial(3)+1 Domain_size(3)]; kz=[-1 0]; kzz=1; 
        Wholevolume(2) = Domain_size(3);
    end
    x0=min(x); x1=max(x);
    y0=min(y); y1=max(y);
    z0=min(z); z1=max(z);    
    
    while 1
        if sub_domain_id>0
            if strcmp(RVEparameters.Growthrelativeto,'% of the previous subvolume')
                current_subdomain_size = [x1-x0+1 y1-y0+1 z1-z0+1];
                dstep = round(RVEparameters.Growthperstep/100 * current_subdomain_size);
            elseif strcmp(RVEparameters.Growthrelativeto,'% of the full volume')
                dstep = round(RVEparameters.Growthperstep/100 * Domain_size);
            end
            x = x+(dstep(1)*kx);
            y = y+(dstep(2)*ky);
            z = z+(dstep(3)*kz);
        end
        x0=min(x); x1=max(x);
        y0=min(y); y1=max(y);
        z0=min(z); z1=max(z);

        if x0<1 || y0<1 || z0<1 || x1>Domain_size(1) || y1>Domain_size(2) || z1>Domain_size(3)
            break % Terminate while loop
        end
        
        sub_domain_group=sub_domain_group+1; % Increment group id
        number_subdomain_group=1; % Number of subdomain
        GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = ( (x1-x0+1)*(y1-y0+1)*(z1-z0+1) ) ^(1/3);
        GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
        GROUP_SUBDOMAIN.id(sub_domain_group).length = abs( (x1-x0+1)*kxx + (y1-y0+1)*kyy + (z1-z0+1)*kzz); % Length;
        GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [(x1-x0+1) (y1-y0+1) (z1-z0+1)]./(z1-z0+1);
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


if strcmp(RVEparameters.type,'H') % 'One subvolume + growing from volume edge (G)'
    d = (RVEparameters.firstuniquevolume_size/100)^(1/2);
    if strcmp(RVEparameters.Constantdirection,'Direction 1')
        Wholevolume(2) = (Domain_size(2) * Domain_size(3))^(1/2);
        x0 = 1; x1 = Domain_size(1);
        y0 = round(Domain_size(2)/2 - d/2*Domain_size(2)); y1 = round(Domain_size(2)/2 + d/2*Domain_size(2));
        z0 = round(Domain_size(3)/2 - d/2*Domain_size(3)); z1 = round(Domain_size(3)/2 + d/2*Domain_size(3));
    elseif strcmp(RVEparameters.Constantdirection,'Direction 2')
        Wholevolume(2) = (Domain_size(1) * Domain_size(3))^(1/2);
        x0 = round(Domain_size(1)/2 - d/2*Domain_size(1)); x1 = round(Domain_size(1)/2 + d/2*Domain_size(1));
        y0 = 1; y1 = Domain_size(2);
        z0 = round(Domain_size(3)/2 - d/2*Domain_size(3)); z1 = round(Domain_size(3)/2 + d/2*Domain_size(3));
    elseif strcmp(RVEparameters.Constantdirection,'Direction 3')
        Wholevolume(2) = (Domain_size(1) * Domain_size(2))^(1/2);
        x0 = round(Domain_size(1)/2 - d/2*Domain_size(1)); x1 = round(Domain_size(1)/2 + d/2*Domain_size(1));
        y0 = round(Domain_size(2)/2 - d/2*Domain_size(2)); y1 = round(Domain_size(2)/2 + d/2*Domain_size(2));
        z0 = 1; z1 = Domain_size(3);
    end
    
    while 1
        if sub_domain_id>0
            g = RVEparameters.Growthperstep/100;
            if strcmp(RVEparameters.Growthrelativeto,'% of the previous subvolume')
                v0 = (x1-x0+1)*(y1-y0+1)*(z1-z0+1); % Current volume
                dv = v0*g;
            elseif strcmp(RVEparameters.Growthrelativeto,'% of the full volume')
                dv = prod(Domain_size)*g;
            end 
            v0 = (x1-x0+1)*(y1-y0+1)*(z1-z0+1); % Current volume
            v1 = v0+dv; % New volume
            d = (v1/prod(Domain_size))^(1/2);
            if strcmp(RVEparameters.Constantdirection,'Direction 1')
                x0 = 1; x1 = Domain_size(1);
                y0 = round(Domain_size(2)/2 - d/2*Domain_size(2)); y1 = round(Domain_size(2)/2 + d/2*Domain_size(2));
                z0 = round(Domain_size(3)/2 - d/2*Domain_size(3)); z1 = round(Domain_size(3)/2 + d/2*Domain_size(3));
            elseif strcmp(RVEparameters.Constantdirection,'Direction 2')
                x0 = round(Domain_size(1)/2 - d/2*Domain_size(1)); x1 = round(Domain_size(1)/2 + d/2*Domain_size(1));
                y0 = 1; y1 = Domain_size(2);
                z0 = round(Domain_size(3)/2 - d/2*Domain_size(3)); z1 = round(Domain_size(3)/2 + d/2*Domain_size(3));
            elseif strcmp(RVEparameters.Constantdirection,'Direction 3')
                x0 = round(Domain_size(1)/2 - d/2*Domain_size(1)); x1 = round(Domain_size(1)/2 + d/2*Domain_size(1));
                y0 = round(Domain_size(2)/2 - d/2*Domain_size(2)); y1 = round(Domain_size(2)/2 + d/2*Domain_size(2));
                z0 = 1; z1 = Domain_size(3);
            end
        end
        if x0<1 || y0<1 || z0<1 || x1>Domain_size(1) || y1>Domain_size(2) || z1>Domain_size(3)
            break % Terminate while loop
        end

        sub_domain_group=sub_domain_group+1; % Increment group id
        number_subdomain_group=1; % Number of subdomain
        GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = ( (x1-x0+1)*(y1-y0+1)*(z1-z0+1) ) ^(1/3);
        if strcmp(RVEparameters.Constantdirection,'Direction 1')
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = ((y1-y0+1)*(z1-z0+1))^(1/2);
        elseif strcmp(RVEparameters.Constantdirection,'Direction 2')
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = ((x1-x0+1)*(z1-z0+1))^(1/2);
        elseif strcmp(RVEparameters.Constantdirection,'Direction 3')
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = ((x1-x0+1)*(y1-y0+1))^(1/2);
        end
        GROUP_SUBDOMAIN.id(sub_domain_group).length = 0; % Length;
        GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [(x1-x0+1) (y1-y0+1) (z1-z0+1)]./(z1-z0+1);
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

end




