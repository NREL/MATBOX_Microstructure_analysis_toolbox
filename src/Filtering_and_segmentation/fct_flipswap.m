function [M,newtype,p] = fct_flipswap(M,p)
sz = size(M);
dimension = length(sz);
newtype = 'same';

new_sz=zeros(1,3);
if dimension==3
    if strcmp(p.action,'Flip axis 1')
        tmp=zeros(sz);
        for k=1:1:sz(1)
            slice=M(k,:,:);
            tmp(sz(1)-k+1,:,:)=slice;
        end
    elseif strcmp(p.action,'Flip axis 2')
        tmp=zeros(sz);
        for k=1:1:sz(2)
            slice=M(:,k,:);
            tmp(:,sz(2)-k+1,:)=slice;
        end
    elseif strcmp(p.action,'Flip axis 3')
        tmp=zeros(sz);
        for k=1:1:sz(3)
            slice=M(:,:,k);
            tmp(:,:,sz(3)-k+1)=slice;
        end
    elseif strcmp(p.action,'Swap axis 1 with axis 2')
        new_sz(1)=sz(2); new_sz(2)=sz(1); new_sz(3)=sz(3);
        tmp=zeros(new_sz);
        for k=1:1:sz(2)
            slice=M(:,k,:);
            tmp(k,:,:)=slice;
        end
    elseif strcmp(p.action,'Swap axis 1 with axis 3')
        new_sz(1)=sz(3); new_sz(2)=sz(2); new_sz(3)=sz(1);
        tmp=zeros(new_sz);
        for k=1:1:sz(3)
            slice=M(:,:,k)';
            tmp(k,:,:)=slice;
        end
    elseif strcmp(p.action,'Swap axis 2 with axis 3')
        new_sz(1)=sz(1); new_sz(2)=sz(3); new_sz(3)=sz(2);
        tmp=zeros(new_sz);
        for k=1:1:sz(3)
            slice=M(:,:,k);
            tmp(:,k,:)=slice;
        end
    end

else
    if strcmp(p.action,'Flip axis 1')
        tmp=zeros(sz);
        for k=1:1:sz(1)
            slice=M(k,:);
            tmp(sz(1)-k+1,:)=slice;
        end
    elseif strcmp(p.action,'Flip axis 2')
        tmp=zeros(sz);
        for k=1:1:sz(2)
            slice=M(:,k);
            tmp(:,sz(2)-k+1)=slice;
        end
    elseif strcmp(p.action,'Swap axis 1 with axis 2')
        new_sz=zeros(1,2); new_sz(1)=sz(2); new_sz(2)=sz(1);
        tmp=zeros(new_sz);
        for k=1:1:sz(2)
            slice=M(:,k);
            tmp(k,:)=slice;
        end
    end

end
M=tmp;

% Convert to uint8 or uint16
[M] = fct_intconvert(M);

end