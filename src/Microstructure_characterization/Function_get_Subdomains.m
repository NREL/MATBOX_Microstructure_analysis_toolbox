function [All_subdomain,GROUP_SUBDOMAIN, Wholevolume] = Function_get_Subdomains(RVEparameters, Domain_size, voxel_size)
% Determine the bounds of subdomains from a domain.

sub_domain_id=0;
sub_domain_group=0;

divisions = RVEparameters.divisions;

voxel_size = voxel_size/1000; % nm -> um

Wholevolume(1) = (prod(Domain_size))^(1/3);
Wholevolume(3) = Domain_size(1)/Domain_size(3); % Aspect ratio
Wholevolume(4) = Domain_size(2)/Domain_size(3);
Wholevolume(5) = Domain_size(3)/Domain_size(3);

if strcmp(RVEparameters.type,'A') % 'Independant subvolumes of same size + keep initial aspect ratio (A)'
    for k=1:1:length(divisions)
        sub_domain_group=sub_domain_group+1;
        number_subdomain_group=divisions(k)^3;
        GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1)/divisions(k) * Domain_size(2)/divisions(k) * Domain_size(3)/divisions(k))^(1/3);
        GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
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
                    All_subdomain(sub_domain_id,6) = x0;
                    All_subdomain(sub_domain_id,7) = x1;
                    All_subdomain(sub_domain_id,8) = y0;
                    All_subdomain(sub_domain_id,9) = y1;
                    All_subdomain(sub_domain_id,10) = z0;
                    All_subdomain(sub_domain_id,11) = z1;
                    All_subdomain(sub_domain_id,12) = divisions(k);
                    
                    tmp_id = [tmp_id sub_domain_id];
                end
            end
        end
        GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
    end
    
    if RVEparameters.subs4
        sub_domain_group=sub_domain_group+1; % Increment group id
        number_subdomain_group=4; % Number of subdomain
        tmp_id = [];        
        if Domain_size(1)<=Domain_size(2) && Domain_size(1)<=Domain_size(3)
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(2)/2 * Domain_size(3)/2 * Domain_size(1))^(1/3);
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
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
                    All_subdomain(sub_domain_id,6) = xx_start(kx);
                    All_subdomain(sub_domain_id,7) = xx_end(kx);
                    All_subdomain(sub_domain_id,8) = yy_start(ky);
                    All_subdomain(sub_domain_id,9) = yy_end(ky);
                    All_subdomain(sub_domain_id,10) = zz_start(kz);
                    All_subdomain(sub_domain_id,11) = zz_end(kz);
                    All_subdomain(sub_domain_id,12) = 0;
                    tmp_id = [tmp_id sub_domain_id];
                end
            end
        end

        GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
    end
    
    if RVEparameters.subs2
        sub_domain_group=sub_domain_group+1; % Increment group id
        number_subdomain_group=2; % Number of subdomain
        tmp_id = [];        
        if Domain_size(1)<=Domain_size(2) && Domain_size(1)<=Domain_size(3)
            if Domain_size(2)>=Domain_size(3)
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(2)/2 * Domain_size(3) * Domain_size(1))^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1) Domain_size(2)/2 Domain_size(3)] ./ Domain_size(3);
                yy_start(1)=1; yy_end(1)=floor(Domain_size(2)/2);
                yy_start(2)=yy_end(1)+1; yy_end(2)=Domain_size(2);
                zz_start(1)=1; zz_end(1)=Domain_size(3);
                zz_start(2)=1; zz_end(2)=Domain_size(3);
            else
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(2) * Domain_size(3)/2 * Domain_size(1))^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
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
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1)/2 Domain_size(2) Domain_size(3)] ./ Domain_size(3);
                xx_start(1)=1; xx_end(1)=floor(Domain_size(1)/2);
                xx_start(2)=xx_end(1)+1; xx_end(2)=Domain_size(1);
                zz_start(1)=1; zz_end(1)=Domain_size(3);
                zz_start(2)=1; zz_end(2)=Domain_size(3);
            else
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1) * Domain_size(3)/2 * Domain_size(2))^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
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
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1)/2 Domain_size(2) Domain_size(3)] ./ Domain_size(3);
                xx_start(1)=1; xx_end(1)=floor(Domain_size(1)/2);
                xx_start(2)=xx_end(1)+1; xx_end(2)=Domain_size(1);
                yy_start(1)=1; yy_end(1)=Domain_size(2);
                yy_start(2)=1; yy_end(2)=Domain_size(2);
            else
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1) * Domain_size(2)/2 * Domain_size(3))^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0;
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
            All_subdomain(sub_domain_id,6) = xx_start(sub2);
            All_subdomain(sub_domain_id,7) = xx_end(sub2);
            All_subdomain(sub_domain_id,8) = yy_start(sub2);
            All_subdomain(sub_domain_id,9) = yy_end(sub2);
            All_subdomain(sub_domain_id,10) = zz_start(sub2);
            All_subdomain(sub_domain_id,11) = zz_end(sub2);
            All_subdomain(sub_domain_id,12) = 0;
            tmp_id = [tmp_id sub_domain_id];
        end
        GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
    end
end


if  strcmp(RVEparameters.type,'B') % 'Independant subvolumes of same size + user-defined aspect ratio (B)'
    
    AR_desired = RVEparameters.Aspectratio / RVEparameters.Aspectratio(3);
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
                    All_subdomain(sub_domain_id,6) = x0;
                    All_subdomain(sub_domain_id,7) = x1;
                    All_subdomain(sub_domain_id,8) = y0;
                    All_subdomain(sub_domain_id,9) = y1;
                    All_subdomain(sub_domain_id,10) = z0;
                    All_subdomain(sub_domain_id,11) = z1;
                    All_subdomain(sub_domain_id,12) = 0;
                    
                    tmp_id = [tmp_id sub_domain_id];
                end
            end
        end
        GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
    end
end


if strcmp(RVEparameters.type,'D') % 'One subvolume + growing from volume center (D)'
    if strcmp(RVEparameters.firstuniquevolume_unit,'% of total volume')
        vinitial = RVEparameters.firstuniquevolume_size * Domain_size / 100;
    elseif strcmp(RVEparameters.firstuniquevolume_unit,'micrometers')
        vinitial = RVEparameters.firstuniquevolume_size * Domain_size / (Domain_size*voxel_size/1000);
    end    
    linitial_AR1 = (prod(vinitial))^(1/3); % Aspect ratio 1 1 1, correct initial size
    linitial = linitial_AR1 * RVEparameters.Aspectratio; % Desired aspect ratio, but wrong initial size
    linitial_AR = linitial / (prod(linitial)^(1/3) / linitial_AR1); % Desired aspect ratio and initial size
    
    x0 = round(Domain_size(1)/2 - linitial_AR(1)/2); x1 = x0+ round(linitial_AR(1)) -1;
    y0 = round(Domain_size(2)/2 - linitial_AR(2)/2); y1 = y0+ round(linitial_AR(2)) -1;
    z0 = round(Domain_size(3)/2 - linitial_AR(3)/2); z1 = z0+ round(linitial_AR(3)) -1;
    
    if strcmp(RVEparameters.Growthrelativeto,'% of total volume') || strcmp(RVEparameters.Growthrelativeto,'micrometers')
        if strcmp(RVEparameters.Growthrelativeto,'% of total volume')
            dstep = (RVEparameters.Growthperstep * Domain_size / 100) /2;
        elseif strcmp(RVEparameters.Growthrelativeto,'micrometers')
            dstep = (RVEparameters.Growthperstep * Domain_size / (Domain_size*voxel_size/1000)) /2;
        end
        dstep_AR1 = (prod(dstep))^(1/3); % Aspect ratio 1 1 1, correct initial size
        dstep = dstep_AR1 * RVEparameters.Aspectratio; % Desired aspect ratio, but wrong initial size
        dstep_AR = round(dstep / (prod(dstep)^(1/3) / dstep_AR1)); % Desired aspect ratio and initial size
        dstep_AR=max(dstep_AR,1); % Minimum step is 1
    end
    
    while 1
        if sub_domain_id>0
            if strcmp(RVEparameters.Growthrelativeto,'% of current subvolume')
                current_subdomain_size = Domain_size; % Initialization
                current_subdomain_size(1) = x1-x0+1;
                current_subdomain_size(2) = y1-y0+1;
                current_subdomain_size(3) = z1-z0+1;
                dstep = (RVEparameters.Growthperstep * current_subdomain_size / 100) /2;
                dstep_AR1 = (prod(dstep))^(1/3); % Aspect ratio 1 1 1, correct initial size
                dstep = dstep_AR1 * RVEparameters.Aspectratio; % Desired aspect ratio, but wrong initial size
                dstep_AR = round(dstep / (prod(dstep)^(1/3) / dstep_AR1)); % Desired aspect ratio and initial size
                dstep_AR=max(dstep_AR,1); % Minimum step is 1
            end
            x0 = x0 - dstep_AR(1); x1 = x1 + dstep_AR(1)-1;
            y0 = y0 - dstep_AR(2); y1 = y1 + dstep_AR(2)-1;
            z0 = z0 - dstep_AR(3); z1 = z1 + dstep_AR(3)-1;
        end
        if x0<1 || y0<1 || z0<1 || x1>Domain_size(1) || y1>Domain_size(2) || z1>Domain_size(3)
            break % Terminate while loop
        end
        
        sub_domain_group=sub_domain_group+1; % Increment group id
        number_subdomain_group=1; % Number of subdomain
        GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = ( (x1-x0+1)*(y1-y0+1)*(z1-z0+1) ) ^(1/3);
        GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = 0; % Length
        GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [(x1-x0+1) (y1-y0+1) (z1-z0+1)]./(z1-z0+1);
        tmp_id = [];
        sub_domain_id = sub_domain_id+1;
        All_subdomain(sub_domain_id,1) = sub_domain_id;
        All_subdomain(sub_domain_id,2) = sub_domain_group;
        All_subdomain(sub_domain_id,3) = number_subdomain_group;
        All_subdomain(sub_domain_id,4) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length;
        All_subdomain(sub_domain_id,5) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length;
        All_subdomain(sub_domain_id,6) = x0;
        All_subdomain(sub_domain_id,7) = x1;
        All_subdomain(sub_domain_id,8) = y0;
        All_subdomain(sub_domain_id,9) = y1;
        All_subdomain(sub_domain_id,10) = z0;
        All_subdomain(sub_domain_id,11) = z1;
        All_subdomain(sub_domain_id,12) = 0;
        tmp_id = [tmp_id sub_domain_id];
        GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
    end    
end


if strcmp(RVEparameters.type,'E') % 'One subvolume + growing from volume edge (E)'
    if strcmp(RVEparameters.firstuniquevolume_unit,'% of total volume')
        dinitial = round(RVEparameters.firstuniquevolume_size * Domain_size / 100);
    elseif strcmp(RVEparameters.firstuniquevolume_unit,'micrometers')
        dinitial = round(RVEparameters.firstuniquevolume_size * Domain_size / (Domain_size*voxel_size/1000));
    end
    
    if strcmp(RVEparameters.Growthrelativeto,'% of total volume')
        dstep = round(RVEparameters.Growthperstep * Domain_size / 100);
    elseif strcmp(RVEparameters.Growthrelativeto,'micrometers')
        dstep = round(RVEparameters.Growthperstep * Domain_size / (Domain_size*voxel_size/1000));
    end
    
    if strcmp(RVEparameters.Growthdirection, 'Direction 1, from x min to x max')
        z=[1 Domain_size(3)]; kz=[0 0]; kzz=0;
        y=[1 Domain_size(2)]; ky=[0 0]; kyy=0;
        x(1)=1;
        x(2)=x(1)+dinitial(1)-1; kx=[0 1]; kxx=1;   
        Wholevolume(2) = Domain_size(1);
        
    elseif strcmp(RVEparameters.Growthdirection,'Direction 1, from x max to x min')
        z=[1 Domain_size(3)]; kz=[0 0]; kzz=0;
        y=[1 Domain_size(2)]; ky=[0 0]; kyy=0;
        x(2)=Domain_size(1);
        x(1)=x(2)-dinitial(1)+1; kx=[-1 0]; kxx=1;     
        Wholevolume(2) = Domain_size(1);
        
    elseif strcmp(RVEparameters.Growthdirection,'Direction 2, from x min to x max')
        z=[1 Domain_size(3)]; kz=[0 0]; kzz=0;
        x=[1 Domain_size(1)]; kx=[0 0]; kxx=0;
        y(1)=1;
        y(2)=y(1)+dinitial(2)-1; ky=[0 1]; kyy=1;       
        Wholevolume(2) = Domain_size(2);
        
    elseif strcmp(RVEparameters.Growthdirection,'Direction 2, from x max to x min')
        z=[1 Domain_size(3)]; kz=[0 0]; kzz=0;
        x=[1 Domain_size(1)]; kx=[0 0]; kxx=0;
        y(2)=Domain_size(2);
        y(1)=y(2)-dinitial(2)+1; ky=[-1 0]; kyy=1;        
        Wholevolume(2) = Domain_size(2);
        
    elseif strcmp(RVEparameters.Growthdirection,'Direction 3, from x min to x max')
        x=[1 Domain_size(1)]; kx=[0 0]; kxx=0;
        y=[1 Domain_size(2)]; ky=[0 0]; kyy=0;
        z(1)=1;
        z(2)=z(1)+dinitial(3)-1; kz=[0 1]; kzz=1; 
        Wholevolume(2) = Domain_size(3);
        
    elseif strcmp(RVEparameters.Growthdirection,'Direction 3, from x max to x min')
        x=[1 Domain_size(1)]; kx=[0 0]; kxx=0;
        y=[1 Domain_size(2)]; ky=[0 0]; kyy=0;
        z(2)=Domain_size(3);
        z(1)=z(2)-dinitial(3)+1; kz=[-1 0]; kzz=1;
        Wholevolume(2) = Domain_size(3);
    end
    x0=min(x); x1=max(x);
    y0=min(y); y1=max(y);
    z0=min(z); z1=max(z);    
    
    while 1
        if sub_domain_id>0
            if strcmp(RVEparameters.Growthrelativeto,'% of current subvolume')
                current_subdomain_size = Domain_size; % Initialization
                current_subdomain_size(1) = x1-x0+1;
                current_subdomain_size(2) = y1-y0+1;
                current_subdomain_size(3) = z1-z0+1;                
                dstep = round(RVEparameters.Growthperstep * current_subdomain_size / 100);
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
        GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = abs( (x1-x0+1)*kxx + (y1-y0+1)*kyy + (z1-z0+1)*kzz); % Length
        GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [(x1-x0+1) (y1-y0+1) (z1-z0+1)]./(z1-z0+1);
        tmp_id = [];
        sub_domain_id = sub_domain_id+1;
        All_subdomain(sub_domain_id,1) = sub_domain_id;
        All_subdomain(sub_domain_id,2) = sub_domain_group;
        All_subdomain(sub_domain_id,3) = number_subdomain_group;
        All_subdomain(sub_domain_id,4) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length;
        All_subdomain(sub_domain_id,5) = GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length;
        All_subdomain(sub_domain_id,6) = x0;
        All_subdomain(sub_domain_id,7) = x1;
        All_subdomain(sub_domain_id,8) = y0;
        All_subdomain(sub_domain_id,9) = y1;
        All_subdomain(sub_domain_id,10) = z0;
        All_subdomain(sub_domain_id,11) = z1;
        All_subdomain(sub_domain_id,12) = 0;
        tmp_id = [tmp_id sub_domain_id];
        GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
    end
end


if strcmp(RVEparameters.type,'C') % 'Independant subvolumes of same size + constant length (C)'
    for k=1:1:length(divisions)
        sub_domain_group=sub_domain_group+1; % Increment group id
        number_subdomain_group=divisions(k)^2; % Number of subdomain
        if strcmp(RVEparameters.Constantdirection,'Direction 1')
            Wholevolume(2) = (Domain_size(2) * Domain_size(3))^(1/2);
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(2)/divisions(k) * Domain_size(3)/divisions(k) * Domain_size(1))^(1/3);
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (Domain_size(2)/divisions(k) * Domain_size(3)/divisions(k))^(1/2);
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
                    All_subdomain(sub_domain_id,6) = x0;
                    All_subdomain(sub_domain_id,7) = x1;
                    All_subdomain(sub_domain_id,8) = y0;
                    All_subdomain(sub_domain_id,9) = y1;
                    All_subdomain(sub_domain_id,10) = z0;
                    All_subdomain(sub_domain_id,11) = z1;
                    All_subdomain(sub_domain_id,12) = divisions(k);
                    
                    tmp_id = [tmp_id sub_domain_id];
                end
            end            
            
        elseif strcmp(RVEparameters.Constantdirection,'Direction 2')
            Wholevolume(2) = (Domain_size(1) * Domain_size(3))^(1/2);
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1)/divisions(k) * Domain_size(3)/divisions(k) * Domain_size(2))^(1/3);
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (Domain_size(1)/divisions(k) * Domain_size(3)/divisions(k))^(1/2);
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
                    All_subdomain(sub_domain_id,6) = x0;
                    All_subdomain(sub_domain_id,7) = x1;
                    All_subdomain(sub_domain_id,8) = y0;
                    All_subdomain(sub_domain_id,9) = y1;
                    All_subdomain(sub_domain_id,10) = z0;
                    All_subdomain(sub_domain_id,11) = z1;
                    All_subdomain(sub_domain_id,12) = divisions(k);
                    
                    tmp_id = [tmp_id sub_domain_id];
                end
            end
            
        elseif strcmp(RVEparameters.Constantdirection,'Direction 3')
            Wholevolume(2) = (Domain_size(1) * Domain_size(2))^(1/2);
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1)/divisions(k) * Domain_size(2)/divisions(k) * Domain_size(3))^(1/3);
            GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (Domain_size(1)/divisions(k) * Domain_size(2)/divisions(k))^(1/2);
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
                    All_subdomain(sub_domain_id,6) = x0;
                    All_subdomain(sub_domain_id,7) = x1;
                    All_subdomain(sub_domain_id,8) = y0;
                    All_subdomain(sub_domain_id,9) = y1;
                    All_subdomain(sub_domain_id,10) = z0;
                    All_subdomain(sub_domain_id,11) = z1;
                    All_subdomain(sub_domain_id,12) = divisions(k);
                    
                    tmp_id = [tmp_id sub_domain_id];
                end
            end
        end
        GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
    end
    if RVEparameters.subs2
        sub_domain_group=sub_domain_group+1; % Increment group id
        number_subdomain_group=2; % Number of subdomain
        tmp_id = [];
        if strcmp(RVEparameters.Constantdirection,'Direction 1')
            if Domain_size(2)>=Domain_size(3)
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(2)/2 * Domain_size(3) * Domain_size(1))^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (Domain_size(2)/2 * Domain_size(3))^(1/2);
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1) Domain_size(2)/2 Domain_size(3)] ./ Domain_size(3);
                yy_start(1)=1; yy_end(1)=floor(Domain_size(2)/2);
                yy_start(2)=yy_end(1)+1; yy_end(2)=Domain_size(2);
                zz_start(1)=1; zz_end(1)=Domain_size(3);
                zz_start(2)=1; zz_end(2)=Domain_size(3);                    
            else
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(2) * Domain_size(3)/2 * Domain_size(1))^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (Domain_size(2) * Domain_size(3)/2)^(1/2);
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
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1)/2 Domain_size(2) Domain_size(3)] ./ Domain_size(3);
                xx_start(1)=1; xx_end(1)=floor(Domain_size(1)/2);
                xx_start(2)=xx_end(1)+1; xx_end(2)=Domain_size(1);
                zz_start(1)=1; zz_end(1)=Domain_size(3);
                zz_start(2)=1; zz_end(2)=Domain_size(3);                
            else
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1) * Domain_size(3)/2 * Domain_size(2))^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (Domain_size(1) * Domain_size(3)/2)^(1/2);
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
                GROUP_SUBDOMAIN.id(sub_domain_group).Aspectratio = [Domain_size(1)/2 Domain_size(2) Domain_size(3)] ./ Domain_size(3);
                xx_start(1)=1; xx_end(1)=floor(Domain_size(1)/2);
                xx_start(2)=xx_end(1)+1; xx_end(2)=Domain_size(1);
                yy_start(1)=1; yy_end(1)=Domain_size(2);
                yy_start(2)=1; yy_end(2)=Domain_size(2);
            else
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_cubic_length = (Domain_size(1) * Domain_size(2)/2 * Domain_size(3))^(1/3);
                GROUP_SUBDOMAIN.id(sub_domain_group).equivalent_square_length = (Domain_size(1) * Domain_size(2)/2)^(1/2);
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
            All_subdomain(sub_domain_id,6) = xx_start(sub2);
            All_subdomain(sub_domain_id,7) = xx_end(sub2);
            All_subdomain(sub_domain_id,8) = yy_start(sub2);
            All_subdomain(sub_domain_id,9) = yy_end(sub2);
            All_subdomain(sub_domain_id,10) = zz_start(sub2);
            All_subdomain(sub_domain_id,11) = zz_end(sub2);
            All_subdomain(sub_domain_id,12) = 0;
            tmp_id = [tmp_id sub_domain_id];
        end
        GROUP_SUBDOMAIN.id(sub_domain_group).domain_id=tmp_id;
    end
end


end




