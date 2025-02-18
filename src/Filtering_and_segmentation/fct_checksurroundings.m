function [Ms,newtype,foo] = fct_checksurroundings(M,p)

sz = size(M);
foo=[];
newtype = 'same';

sz = size(M);
dimension = length(sz);

Ms = M;

idx = find(M==p.label);
if ~isempty(idx)
    if dimension == 2
        [IX,IY] = ind2sub(sz,idx);
    else
        [IX,IY,IZ] = ind2sub(sz,idx);
    end

    for k=1:1:length(idx)
        xmin = max([1 IX(k)-p.windowradius]); xmax = min([sz(1) IX(k)+p.windowradius]);
        ymin = max([1 IY(k)-p.windowradius]); ymax = min([sz(2) IY(k)+p.windowradius]);
        zmin = 1; zmax = 1;
        zlabel = 1;
        if dimension == 3
            zmin = max([1 IZ(k)-p.windowradius]); zmax = min([sz(3) IZ(k)+p.windowradius]);
            if zmin > IZ(k)-p.windowradius
                zlabel = IZ(k);
            else
                zlabel = p.windowradius+1;
            end
        end

        if xmin > IX(k)-p.windowradius
            xlabel = IX(k);
        else
            xlabel = p.windowradius+1;
        end
        if ymin > IY(k)-p.windowradius
            ylabel = IY(k);
        else
            ylabel = p.windowradius+1;
        end
           
        sub1 = M(xmin:xmax, ymin:ymax, zmin:zmax);
        sub1 = sub1+1;

        sub2 = zeros(size(sub1)+2);
        if dimension == 2
            sub2(2:end-1,2:end-1)=sub1;
        else
            sub2(2:end-1,2:end-1,2:end-1)=sub1;
        end

        unis = unique(sub2);
        if max(unis==p.surround_label+1) % Possible
            % Select cluster
            BW = sub2;
            BW(BW~=p.label+1)=0;
            if dimension == 2
                C = bwlabel(BW,4);
            else
                C = bwlabeln(BW,6);
            end
            unis_C = unique(C);
            if length(unis_C)>2
                % HERE
                Clabel = C(xlabel+1,ylabel+1,zlabel+1);
                C(C==0)=Clabel;
                sub2(C~=Clabel)=0;
            end
            
            % Calculate interface will other phases
            otherphases = unis;
            otherphases(otherphases==p.label+1)=[];
            nothers = length(otherphases);
            surfaces = zeros(nothers,1);
            idlabel = find(sub2==p.label+1);
            for kother = 1:1:nothers
                BW = zeros(size(sub2));
                BW(idlabel)=1;
                BW(sub2==otherphases(kother))=10;
                [surfaces(kother)] = Function_Specificinterface_direct_Algorithm_noPos(BW);
                if p.surround_label+1 == otherphases(kother)
                    s_check = surfaces(kother);
                end
            end
            s_tot = sum(surfaces);
            ratio = s_check/s_tot;

            if ratio >= p.thresholdratio
                Ms(idx(k)) = p.assignto;
            end          
            
        end

    end
end

[Ms] = fct_intconvert(Ms);

end