function [M,newtype,diameter] = fct_localdilation_spheresizefilter(M,p)

sz = size(M);
dimension = length(sz);
newtype = 'same';

if p.perslice && dimension==3
    if strcmp(p.alongaxe,'along axe 3')
        diameter = zeros(sz);
        for z=1:1:sz(3)
            [M(:,:,z),diameter(:,:,z)] = fct_localdilation_spheresizefilter_algo(M(:,:,z),p);
        end

    elseif strcmp(p.alongaxe,'along axe 1')
        % swap axis 1 with 3
        new_sz(1)=sz(3); new_sz(2)=sz(2); new_sz(3)=sz(1);
        tmp=zeros(new_sz);
        for k=1:1:sz(3)
            slice=M(:,:,k)';
            tmp(k,:,:)=slice;
        end
        M=tmp;

        sz = size(M);
        diameter = zeros(sz);
        for z=1:1:sz(3)
            [M(:,:,z),diameter(:,:,z)] = fct_localdilation_spheresizefilter_algo(M(:,:,z),p);
        end

        % swap axis 1 with 3
        new_sz(1)=sz(3); new_sz(2)=sz(2); new_sz(3)=sz(1);
        tmp=zeros(new_sz);
        for k=1:1:sz(3)
            slice=M(:,:,k)';
            tmp(k,:,:)=slice;
        end
        M=tmp;
        tmp=zeros(new_sz);
        for k=1:1:sz(3)
            slice=diameter(:,:,k)';
            tmp(k,:,:)=slice;
        end      
        diameter=tmp;  

    elseif strcmp(p.alongaxe,'along axe 2')
        % swap axis 2 with 3
        new_sz(1)=sz(1); new_sz(2)=sz(3); new_sz(3)=sz(2);
        tmp=zeros(new_sz);
        for k=1:1:sz(3)
            slice=M(:,:,k);
            tmp(:,k,:)=slice;
        end
        M=tmp;

        sz = size(M);
        diameter = zeros(sz);
        for z=1:1:sz(3)
            [M(:,:,z),diameter(:,:,z)] = fct_localdilation_spheresizefilter_algo(M(:,:,z),p);
        end

        % swap axis 2 with 3
        new_sz(1)=sz(1); new_sz(2)=sz(3); new_sz(3)=sz(2);
        tmp=zeros(new_sz);
        for k=1:1:sz(3)
            slice=M(:,:,k);
            tmp(:,k,:)=slice;
        end
        M=tmp;
        tmp=zeros(new_sz);
        for k=1:1:sz(3)
            slice=diameter(:,:,k);
            tmp(:,k,:)=slice;
        end
        diameter=tmp;

    end
else
    [M,diameter] = fct_localdilation_spheresizefilter_algo(M,p);
end

% Convert to uint8 or uint16
[M] = fct_intconvert(M);

end