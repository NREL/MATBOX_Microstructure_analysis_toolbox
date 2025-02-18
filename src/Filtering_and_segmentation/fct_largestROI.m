function [M,newtype,res] = fct_largestROI(M,p)
sz = size(M);
dimension = length(sz);
newtype = 'same';

BW = zeros(sz);
if strcmp(p.selectbackground,'From label')
    BW(M~=p.background_label)=1;
else
    BW(p.background_volume~=1)=1;
end

%Add downscaling for dim3 only, increase range

if strcmp(p.ops,'Cropping')
    [x0, x1, y0, y1, z0, z1, ~] = fct_findlargestrectangle3D(BW,p.buffers);
    if dimension==2
        res.crop = [x0 x1 y0 y1];
        res.str_p = ['Axe 1: ' num2str(x0,'%i') '-' num2str(x1,'%i') ', axe 2: ' num2str(y0,'%i') '-' num2str(y1,'%i')];
    else
        res.crop = [x0 x1 y0 y1 z0 z1];
        res.str_p = ['Axe 1: ' num2str(x0,'%i') '-' num2str(x1,'%i') ', axe 2: ' num2str(y0,'%i') '-' num2str(y1,'%i') ', axe 3: ' num2str(z0,'%i') '-' num2str(z1,'%i')];
    end
    if p.apply
        if dimension==2
            M = M(x0:x1,y0:y1);
        else
            M = M(x0:x1,y0:y1,z0:z1);
        end
    end

elseif strcmp(p.ops,'Cropping and rotations')

    n_rot = 5; % Must be odd number
    n_refine = 10;

    candidates = zeros(n_refine,4);

    if dimension == 3
        repeating_background_z = true;
        repeating_background_x = false;
        repeating_background_y = false;

        for z=1:1:sz(3)-1
            if sum(sum( BW(:,:,z) == BW(:,:,z+1) ))~=sz(1)*sz(2)
                repeating_background_z = false;
                break
            end
        end
        if ~repeating_background_z
            repeating_background_x = true;
            for x=1:1:sz(1)-1
                if sum(sum( BW(x,:,:) == BW(x+1,:,:) ))~=sz(2)*sz(3)
                    repeating_background_x = false;
                    break
                end
            end
        end
        if ~repeating_background_z && ~repeating_background_x
            repeating_background_y = true;
            for y=1:1:sz(2)-1
                if sum(sum( BW(:,y,:) == BW(:,y+1,:) ))~=sz(1)*sz(3)
                    repeating_background_y = false;
                    break
                end
            end
        end
    else
        p.rotmax_x = 0;
        p.rotmax_y = 0;
    end

    for k_refine = 1:1:n_refine
        if k_refine==1
            rotxs = unique(linspace(-p.rotmax_x,p.rotmax_x,n_rot));
            rotys = unique(linspace(-p.rotmax_y,p.rotmax_y,n_rot));
            rotzs = unique(linspace(-p.rotmax_z,p.rotmax_z,n_rot));
        else
            rotxs = unique(linspace(-p.rotmax_x/k_refine,p.rotmax_x/k_refine,n_rot))+best_rx;
            rotys = unique(linspace(-p.rotmax_y/k_refine,p.rotmax_y/k_refine,n_rot))+best_ry;
            rotzs = unique(linspace(-p.rotmax_z/k_refine,p.rotmax_z/k_refine,n_rot))+best_rz;
        end

        if dimension == 3
            if repeating_background_x
                rotys = 0; rotzs= 0;
            elseif repeating_background_y
                rotxs = 0; rotzs= 0;
            elseif repeating_background_z
                rotxs = 0; rotys= 0;
            end
        end

        solutions = [];
        pp.addition = false;
        for kx = 1:1:length(rotxs)
            rx = rotxs(kx);
            for ky = 1:1:length(rotys)
                ry = rotys(ky);
                for kz = 1:1:length(rotzs)
                    rz = rotzs(kz);

                    BWrotated = BW;
                    if rx~=0
                        pp.axis = 'Axis 1';
                        pp.angle = rx;
                        [BWrotated,~,~] = fct_rotation(BWrotated,pp);
                    end
                    if ry~=0
                        pp.axis = 'Axis 2';
                        pp.angle = ry;
                        [BWrotated,~,~] = fct_rotation(BWrotated,pp);
                    end
                    if rz~=0
                        pp.axis = 'Axis 3';
                        pp.angle = rz;
                        [BWrotated,~,~] = fct_rotation(BWrotated,pp);
                    end

                    idx = find(BWrotated);
                    if dimension==3
                        [IX, IY, IZ] = ind2sub(size(BWrotated),idx);
                        x0 = min(IX); x1=max(IX);
                        y0 = min(IY); y1=max(IY);
                        z0 = min(IZ); z1=max(IZ);
                        BWrotated = BWrotated(x0:x1,y0:y1,z0:z1); % Otherwise volume growths unnecessarily
                    else
                        [IX, IY] = ind2sub(size(BWrotated),idx);
                        x0 = min(IX); x1=max(IX);
                        y0 = min(IY); y1=max(IY);
                        BWrotated = BWrotated(x0:x1,y0:y1); % Otherwise volume growths unnecessarily
                    end

                    [~, ~, ~, ~, ~, ~, maxvolume] = fct_findlargestrectangle3D(double(BWrotated),p.buffers);
                    solutions = [solutions; rx ry rz maxvolume];
                end
            end
        end
        solutions = sortrows(solutions,4,'descend');
        best_rx = solutions(1,1);
        best_ry = solutions(1,2);
        best_rz = solutions(1,3);
        candidates(k_refine,:) = solutions(1,:);
    end

    candidates = sortrows(candidates,4,'descend');

    res.str_p = 'In this order:';

    % We have to redo because we cropped: BWrotated = BWrotated(x0:x1,y0:y1,z0:z1), and thus we do not have the x0,x1,y0,y1,z0,z1
    rx = candidates(1,1); ry = candidates(1,2);  rz = candidates(1,3);
    if rx~=0
        pp.axis = 'Axis 1';
        pp.angle = rx;
        [BW,~,~] = fct_rotation(BW,pp);
        if p.apply
            [M,~,~] = fct_rotation(M,pp);
        end
        res.str_p = [res.str_p ' Axis of rotation: 1, angle:' num2str(pp.angle,'%1.2f')];
    end
    if ry~=0
        pp.axis = 'Axis 2';
        pp.angle = ry;
        [BW,~,~] = fct_rotation(BW,pp);
        if p.apply
            [M,~,~] = fct_rotation(M,pp);
        end
        res.str_p = [res.str_p ' Axis of rotation: 2, angle:' num2str(pp.angle,'%1.2f')];
    end
    if rz~=0
        pp.axis = 'Axis 3';
        pp.angle = rz;
        [BW,~,~] = fct_rotation(BW,pp);
        if p.apply
            [M,~,~] = fct_rotation(M,pp);
        end
        res.str_p = [res.str_p ' Axis of rotation: 3, angle:' num2str(pp.angle,'%1.2f')];
    end

    [x0, x1, y0, y1, z0, z1, ~] = fct_findlargestrectangle3D(double(BW),p.buffers);

    if dimension==2
        res.str_p = [res.str_p ', Cropping: Axe 1: ' num2str(x0,'%i') '-' num2str(x1,'%i') ', axe 2: ' num2str(y0,'%i') '-' num2str(y1,'%i')];
    else
        res.str_p = [res.str_p ', Cropping: Axe 1: ' num2str(x0,'%i') '-' num2str(x1,'%i') ', axe 2: ' num2str(y0,'%i') '-' num2str(y1,'%i') ', axe 3: ' num2str(z0,'%i') '-' num2str(z1,'%i')];
    end
    if p.apply
        if dimension==2
            M = M(x0:x1,y0:y1);
        else
            M = M(x0:x1,y0:y1,z0:z1);
        end
    end

end

end