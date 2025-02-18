function [x0, x1, y0, y1, z0, z1, maxvolume] = fct_findlargestrectangle3D(M,background_label,buffer)
% Find largest rectangle not withinn the background

sz = size(M);
dimension = length(sz);
if dimension == 2 % 2D: histogram-based algorithm
    BW = zeros(sz);
    BW( M~=background_label ) = 1;
    [x0, x1, y0, y1, ~] = fct_findlargestrectangle2D(BW); 
    if x0+buffer < x1 - buffer
        x0 = x0 + buffer;
        x1 = x1 - buffer;
    end
    if y0+buffer < y1 - buffer
        y0 = y0 + buffer;
        y1 = y1 - buffer;
    end
    maxvolume = (x1-x0+1)*(y1-y0+1);
    z0 = 1;
    z1 = 1;

else % 3D: guided ray-tracing approach. May not find largest rectangle each time.
    BW = zeros(sz);
    BW( M~=background_label ) = 1;
    dmap = bwdist(~BW,'Chessboard');
    idx = find(BW);
    [IX,IY,IZ] = ind2sub(sz,idx);

    figure; imagesc(BW(:,:,round(sz(3)/2))); axis equal; axis tight; colormap gray;
    figure; imagesc(dmap(:,:,round(sz(3)/2))); axis equal; axis tight; colormap turbo;

    choices = [IX IY IZ dmap(idx)];
    choices = sortrows(choices,4,'descend');
    n = length(idx);

    current_largest_rectangle = 0;
    for k=1:1:n
        x = choices(k,1);
        y = choices(k,2);
        z = choices(k,3);
        % % x,y,z probe lines
        xs = BW(:,y,z);
        ys = BW(x,:,z);
        zs = BW(x,y,:);
        zs = reshape(zs,[numel(zs),1]);
        % Find extremities
        id0x = find(xs==0);
        id0y = find(ys==0);
        id0z = find(zs==0);
        [x0,x1] = find_extremities(id0x,x,sz(1))
        [y0,y1] = find_extremities(id0y,y,sz(2))
        [z0,z1] = find_extremities(id0z,z,sz(3))

        xl0 = x0; xl1 = x1;
        yl0 = y0; yl1 = y1;
        zl0 = z0; zl1 = z1;        

        keyboard

        % Largest possible rectangle considering only x,y,z probe lines
        largest_possible_rectangle = (z1-z0+1) * (y1-y0+1) * (x1-x0+1);
        if largest_possible_rectangle<=current_largest_rectangle % Move to next point
            continue
        else % It is possible we can find a larger rectange
            nangle = 1;
            div = 1/2 * 1/(nangle+1);

            % z slice
            for sens=1:1:2
                for angle = -nangle:1:nangle
                    den = 1/div * abs(angle);

                    if angle == 0
                        xx = [x0 x1];
                        yy = [y0 y1];
                    elseif angle <0
                        
                        yy = [y0+(y1-y0)/den y1-(y1-y0)/den];
                        xx = [x0 x1];
                    else
                        xx = [x0+(x1-x0)/den x1-(x1-x0)/den];
                        yy = [y0 y1];
                    end
                    if sens==2
                        yy = flip(yy);
                    end
                    xx = round(xx); yy = round(yy);

                    [xs,ys] = Find_allindex_2Dline(xx(1),xx(2),yy(1),yy(2));
                    idx = sub2ind(sz,xs,ys,z*ones(1,length(xs)));
                    vals = BW(idx);
                    ids = find(vals==0);
                    if ~isempty(ids)
                        idx = find(xs==x);
                        [xa,xb] = find_extremities(ids,idx(1),sz(1));
                        xa = xs(xa); xb = xs(xb);
                        idy = find(ys==y);
                        [ya,yb] = find_extremities(ids,idy(1),sz(2));
                        ya = ys(ya); yb = ys(yb);
                        xl0 = max([xl0 xa]); xl1 = min([xl1 xb]);
                        yl0 = max([yl0 ya]); yl1 = min([yl1 yb]);
                    end

                    BW2=BW(:,:,1);
                    for kk=1:1:length(xs)
                        BW2(xs(kk),ys(kk))=2;
                    end
                    figure; imagesc(BW2); axis equal; axis tight; colormap gray;
                    keyboard

                    [xl0 xl1]
                    [yl0 yl1]

                end
            end




            % 1st diagonale
            [xs,ys] = Find_allindex_2Dline(x0,x1,y0,y1);
            idx = sub2ind(sz,xs,ys,z*ones(1,length(ys)));
            vals = BW(idx);
            ids = find(vals==0);
            idx = find(xs==x);
            [xa,xb] = find_extremities(ids,idx,sz(1));
            xa = xs(xa); xb = xs(xb);
            idy = find(ys==y);
            [ya,yb] = find_extremities(ids,idy,sz(2));
            ya = ys(ya); yb = ys(yb);
            x0 = max([x0 xa]); x1 = min([x1 xb]);
            y0 = max([y0 ya]); y1 = min([y1 yb]);            

            % 2nd diagonale
            [xs,ys] = Find_allindex_2Dline(x0,x1,y1,y0);
            idx = sub2ind(sz,xs,ys,z*ones(1,length(ys)));
            vals = BW(idx);
            ids = find(vals==0);
            idx = find(xs==x);
            [xa,xb] = find_extremities(ids,idx,sz(1));
            xa = xs(xa); xb = xs(xb);
            idy = find(ys==y);
            [ya,yb] = find_extremities(ids,idy,sz(2));
            ya = ys(ya); yb = ys(yb);
            x0 = max([x0 xa]); x1 = min([x1 xb]);
            y0 = max([y0 ya]); y1 = min([y1 yb]);               


            BW2=BW(:,:,1);
            for kk=1:1:length(xs)
                BW2(xs(kk),ys(kk))=2;
            end
            figure; imagesc(BW2); axis equal; axis tight; colormap gray;




            % % x slice
            % 1st diagonale
            [ys,zs] = Find_allindex_2Dline(y0,y1,z0,z1);
            idx = sub2ind(sz,x*ones(1,length(ys)),ys,zs);
            vals = BW(idx);
            ids = find(vals==0);
            idy = find(ys==y);
            [ya,yb] = find_extremities(ids,idy,sz(2))
            idz = find(zs==z);
            [za,zb] = find_extremities(ids,idz,sz(3))
            y0 = max([y0 ya]); y1 = min([y1 yb]);
            z0 = max([z0 za]); z1 = min([z1 zb]);
            % 2nd diagonale
            [ys,zs] = Find_allindex_2Dline(y0,y1,z1,z0);
            idx = sub2ind(sz,x*ones(1,length(ys)),ys,zs);
            vals = BW(idx);
            ids = find(vals==0);
            [ya,yb] = find_extremities(ids,idy,sz(2))
            idz = find(zs==z);
            [za,zb] = find_extremities(ids,idz,sz(3))
            y0 = max([y0 ya]); y1 = min([y1 yb]);
            z0 = max([z0 za]); z1 = min([z1 zb]);

            x
            [y0 y1]
            [z0 z1]

            llll


        end

                

        % % Largest possible rectangle considering only x,y,z probe lines
        % largest_possible_rectangle = (z_max-z_min+1) * (y_max-y_min+1) * (x_max-x_min+1);
        % if largest_possible_rectangle<=current_largest_rectangle % Move to next point
        %     continue 
        % else % It is possible we can find a larger rectange
        %     for x0 = x_min:1:x_max-1
        %         for x1 = x0+1:1:x_max
        %             for y0 = y_min:1:y_max-1
        %                 for y1 = y0+1:1:y_max
        %                     for z0 = z_min:1:z_max-1
        %                         for z1 = z0+1:1:z_max
        %                             sub = BW(x0:x1,y0:y1,z0:z1);
        %                             n0 = sum(sum(sum(sub==0)));
        %                             if n0==0
        %                                 v = (z1-z0+1) * (y1-y0+1) * (x1-x0+1);
        %                                 if v>current_largest_rectangle
        %                                     current_largest_rectangle = v;
        %                                     current_x0 = x0;
        %                                     current_x1 = x1;
        %                                     current_y0 = y0;
        %                                     current_y1 = y1;
        %                                     current_z0 = z0;
        %                                     current_z1 = z1;
        %                                 end
        %                             end                               
        %                         end
        %                     end
        %                 end
        %             end
        %         end
        %     end
        % 
        % end

    end

end