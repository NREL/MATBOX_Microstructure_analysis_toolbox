function [x0, x1, y0, y1, z0, z1, maxvolume] = fct_findlargestrectangle3D(BW,buffers)
% Find largest rectangle not withinn the background

sz = size(BW);
dimension = length(sz);
if dimension == 2 % 2D: histogram-based algorithm
    [x0, x1, y0, y1, ~] = fct_findlargestrectangle2D(BW);
    if x0+buffers(1) < x1 - buffers(1)
        x0 = x0 + buffers(1);
        x1 = x1 - buffers(1);
    end
    if y0+buffers(2) < y1 - buffers(2)
        y0 = y0 + buffers(2);
        y1 = y1 - buffers(2);
    end
    maxvolume = (x1-x0+1)*(y1-y0+1);
    z0 = 1;
    z1 = 1;

else % 3D

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

    if repeating_background_z % Same with 2D case
        [x0, x1, y0, y1, ~] = fct_findlargestrectangle2D(BW(:,:,1));
        if x0+buffers(1) < x1 - buffers(1)
            x0 = x0 + buffers(1);
            x1 = x1 - buffers(1);
        end
        if y0+buffers(2) < y1 - buffers(2)
            y0 = y0 + buffers(2);
            y1 = y1 - buffers(2);
        end
        z0 = 1; z1 = sz(3);
        if z0+buffers(3) < z1 - buffers(3)
            z0 = z0 + buffers(3);
            z1 = z1 - buffers(3);
        end
        maxvolume = (x1-x0+1)*(y1-y0+1)*(z1-z0+1);
    
    elseif repeating_background_y % Same with 2D case
        sl = squeeze(BW(:,1,:));
        [x0, x1, z0, z1, ~] = fct_findlargestrectangle2D(sl);
        if x0+buffers(1) < x1 - buffers(1)
            x0 = x0 + buffers(1);
            x1 = x1 - buffers(1);
        end
        if z0+buffers(3) < z1 - buffers(3)
            z0 = z0 + buffers(3);
            z1 = z1 - buffers(3);
        end
        y0 = 1; y1 = sz(2);
        if y0+buffers(2) < y1 - buffers(2)
            y0 = y0 + buffers(2);
            y1 = y1 - buffers(2);
        end
        maxvolume = (x1-x0+1)*(y1-y0+1)*(z1-z0+1);

    elseif repeating_background_x % Same with 2D case
        sl = squeeze(BW(1,:,:));
        [y0, y1, z0, z1, ~] = fct_findlargestrectangle2D(sl);
        if y0+buffers(2) < y1 - buffers(2)
            y0 = y0 + buffers(2);
            y1 = y1 - buffers(2);
        end
        if z0+buffers(3) < z1 - buffers(3)
            z0 = z0 + buffers(3);
            z1 = z1 - buffers(3);
        end
        x0 = 1; x1 = sz(1);
        if x0+buffers(1) < x1 - buffers(1)
            x0 = x0 + buffers(1);
            x1 = x1 - buffers(1);
        end
        maxvolume = (x1-x0+1)*(y1-y0+1)*(z1-z0+1);        

    else % Probability approach. May not find largest rectangle each time.
        autosdownscale = true;
        max_numel = 150e3;
        % Probability step
        n_refinement = 5; % 10
        n_samples_equi = 1000; % 5000
        n_samples_prob = 500; % 1000
        % Refinment step
        n_iterate = 4;
        range = 2;

        if numel(BW)<2*max_numel
            autosdownscale = false; % Not worth it
        end
        if autosdownscale
            sz_original = sz;
            BW_original = BW;
            pp.scaling_factor = (numel(BW)/max_numel)^(1/3);
            pp.label_or_greylevel = 'Label';
            pp.background = 0;
            [BW] = function_scaling(BW,pp);
            sz = size(BW);
        end

        largest_restangle = 0;
        allborders = [];
        weights = [];
        %allres = [];
        for k_refinemnt = 1:1:n_refinement
            if k_refinemnt==1 % Equiprobability
                xx = randi(sz(1),n_samples_equi,2);
                xx = sort(xx,2);
                yy = randi(sz(2),n_samples_equi,2);
                yy = sort(yy,2);
                zz = randi(sz(3),n_samples_equi,2);
                zz = sort(zz,2);
                n_samples_current = n_samples_equi;
            else % Calculate probability based on all previous attemps and samples within this distribution
                [x0s] = sample_from_distributions(sz(1),allborders(:,1),weights,n_samples_prob);
                [x1s] = sample_from_distributions(sz(1),allborders(:,2),weights,n_samples_prob);
                [y0s] = sample_from_distributions(sz(2),allborders(:,3),weights,n_samples_prob);
                [y1s] = sample_from_distributions(sz(2),allborders(:,4),weights,n_samples_prob);
                [z0s] = sample_from_distributions(sz(3),allborders(:,5),weights,n_samples_prob);
                [z1s] = sample_from_distributions(sz(3),allborders(:,6),weights,n_samples_prob);
                xx = round([x0s' x1s']);
                yy = round([y0s' y1s']);
                zz = round([z0s' z1s']);
                xx(xx(:,2)<xx(:,1),:)=[];
                yy(yy(:,2)<yy(:,1),:)=[];
                zz(zz(:,2)<zz(:,1),:)=[];
                [nx,~]=size(xx);
                [ny,~]=size(yy);
                [nz,~]=size(zz);
                n_samples_current = min([nx ny nz]);
            end

            for k=1:1:n_samples_current
                x0 = xx(k,1); x1 = xx(k,2);
                y0 = yy(k,1); y1 = yy(k,2);
                z0 = zz(k,1); z1 = zz(k,2);
                BWtmp = BW(x0:x1,y0:y1,z0:z1);
                if ~any(~BWtmp,'all')
                    current_rectangle = (z1-z0+1) * (y1-y0+1) * (x1-x0+1);
                    allborders = [allborders; x0 x1 y0 y1 z0 z1];
                    weights = [weights; current_rectangle];
                    if current_rectangle>largest_restangle
                        largest_restangle = current_rectangle;
                        best = [x0 x1 y0 y1 z0 z1];
                        %allborders = [allborders; x0 x1 y0 y1 z0 z1 current_rectangle];
                    end
                end
            end

            %allres(k_refinemnt) = largest_restangle;
        end
        %figure; plot([1:1:n_refinement],allres)

        % Iterate along the edges. Very costly, so range should be very low.
        % This step will refine the solution from the solution found with the probabilistic approach
        % that is we are already near the solution, so range does not have to be high anyway
        for k_iterate = 1:1:n_iterate
            x0s = [best(1)-range:1: best(1)+range];
            x0s(x0s<1)=[]; x0s(x0s>sz(1))=[];
            x1s = [best(2)-range:1: best(2)+range];
            x1s(x1s<1)=[]; x1s(x1s>sz(1))=[];

            y0s = [best(3)-range:1: best(3)+range];
            y0s(y0s<1)=[]; y0s(y0s>sz(2))=[];
            y1s = [best(4)-range:1: best(4)+range];
            y1s(y1s<1)=[]; y1s(y1s>sz(2))=[];

            z0s = [best(5)-range:1: best(5)+range];
            z0s(z0s<1)=[]; z0s(z0s>sz(3))=[];
            z1s = [best(6)-range:1: best(6)+range];
            z1s(z1s<1)=[]; z1s(z1s>sz(3))=[];

            for kx0 = 1:1:length(x0s)
                x0 = x0s(kx0);
                for kx1 = 1:1:length(x1s)
                    x1 = x1s(kx1);
                    for ky0 = 1:1:length(y0s)
                        y0 = y0s(ky0);
                        for ky1 = 1:1:length(y1s)
                            y1 = y1s(ky1);
                            for kz0 = 1:1:length(z0s)
                                z0 = z0s(kz0);
                                for kz1 = 1:1:length(z1s)
                                    z1 = z1s(kz1);

                                    BWtmp = BW(x0:x1,y0:y1,z0:z1);
                                    if ~any(~BWtmp,'all')
                                        current_rectangle = (z1-z0+1) * (y1-y0+1) * (x1-x0+1);
                                        if current_rectangle>largest_restangle
                                            largest_restangle = current_rectangle;
                                            best_new = [x0 x1 y0 y1 z0 z1];
                                        end
                                    end

                                end
                            end
                        end
                    end
                end
            end

            if sum(best_new==best)==6
                break
            else
                best = best_new;
            end
        end
        x0 = best_new(1); x1 = best_new(2);
        y0 = best_new(3); y1 = best_new(4);
        z0 = best_new(5); z1 = best_new(6);

        if autosdownscale
            if x0==1
                x0=1;
            else
                x0 = round(x0*pp.scaling_factor) + 1;
            end
            if x1==sz(1)
                x1=sz_original(1);
            else
                x1 = round((x1-1)*pp.scaling_factor);
            end
            if y0==1
                y0=1;
            else
                y0 = round(y0*pp.scaling_factor) + 1;
            end
            if y1==sz(2)
                y1=sz_original(2);
            else
                y1 = round((y1-1)*pp.scaling_factor);
            end         
            if z0==1
                z0=1;
            else
                z0 = round(z0*pp.scaling_factor) + 1;
            end
            if z1==sz(3)
                z1=sz_original(3);
            else
                z1 = round((z1-1)*pp.scaling_factor);
            end             
            x0 = max([1 x0]); y0 = max([1 y0]); z0 = max([1 z0]);
            x1 = min([sz_original(1) x1]); y1 = min([sz_original(2) y1]); z1 = min([sz_original(3) z1]);

            % best = [x0 x1 y0 y1 z0 z1]
            % largest_restangle = 0;
            % 
            % range = 0;
            % while largest_restangle == 0
            %     range = range + 1
            %     keyboard
            % 
            %     x0s = [best(1)-range:1: best(1)+range];
            %     x0s(x0s<1)=[]; x0s(x0s>sz_original(1))=[];
            %     x1s = [best(2)-range:1: best(2)+range];
            %     x1s(x1s<1)=[]; x1s(x1s>sz_original(1))=[];
            % 
            %     y0s = [best(3)-range:1: best(3)+range];
            %     y0s(y0s<1)=[]; y0s(y0s>sz_original(2))=[];
            %     y1s = [best(4)-range:1: best(4)+range];
            %     y1s(y1s<1)=[]; y1s(y1s>sz_original(2))=[];
            % 
            %     z0s = [best(5)-range:1: best(5)+range];
            %     z0s(z0s<1)=[]; z0s(z0s>sz_original(3))=[];
            %     z1s = [best(6)-range:1: best(6)+range];
            %     z1s(z1s<1)=[]; z1s(z1s>sz_original(3))=[];
            % 
            %     for kx0 = 1:1:length(x0s)
            %         x0 = x0s(kx0)
            %         for kx1 = 1:1:length(x1s)
            %             x1 = x1s(kx1);
            %             for ky0 = 1:1:length(y0s)
            %                 y0 = y0s(ky0);
            %                 for ky1 = 1:1:length(y1s)
            %                     y1 = y1s(ky1);
            %                     for kz0 = 1:1:length(z0s)
            %                         z0 = z0s(kz0);
            %                         for kz1 = 1:1:length(z1s)
            %                             z1 = z1s(kz1);
            % 
            %                             BWtmp = BW_original(x0:x1,y0:y1,z0:z1);
            %                             if ~any(~BWtmp,'all')
            %                                 current_rectangle = (z1-z0+1) * (y1-y0+1) * (x1-x0+1);
            %                                 if current_rectangle>largest_restangle
            %                                     largest_restangle = current_rectangle;
            %                                     best_new = [x0 x1 y0 y1 z0 z1];
            %                                 end
            %                             end                                     
            % 
            %                         end
            %                     end
            %                 end
            %             end
            %         end
            %     end
            % 
            %     toc
            % 
            % end
            % x0 = best_new(1); x1 = best_new(2);
            % y0 = best_new(3); y1 = best_new(4);
            % z0 = best_new(5); z1 = best_new(6);
        end

        if x0+buffers(1) < x1 - buffers(1)
            x0 = x0 + buffers(1);
            x1 = x1 - buffers(1);
        end
        if y0+buffers(2) < y1 - buffers(2)
            y0 = y0 + buffers(2);
            y1 = y1 - buffers(2);
        end
        if z0+buffers(3) < z1 - buffers(3)
            z0 = z0 + buffers(3);
            z1 = z1 - buffers(3);
        end

        maxvolume = (x1-x0+1)*(y1-y0+1)*(z1-z0+1);

    end

end

    function [res] = sample_from_distributions(n,vals,weights,n_samples)
        x_ = 1:1:n;
        pdf = zeros(1,n);
        tot_w = sum(weights);
        for k_ = 1:1:n
            idx_ = find(vals==x_(k_));
            pdf(k_)= sum( weights(idx_) ) / tot_w;
        end
        % Cumulative function
        cdf = zeros(1,n);
        cdf(end) = pdf(end);
        for k_=n-1:-1:1
            cdf(k_) = cdf(k_+1) + pdf(k_) +1e-9; % 1e-9 to avoid identic value that would prevent interp1
        end
        % Sample from distribution
        res = interp1(cdf,x_,rand(1,n_samples),'nearest','extrap');
    end

end
