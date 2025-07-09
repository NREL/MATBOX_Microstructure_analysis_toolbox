function [M,newtype,foo] = fct_remvovegapinnumbering(M,p)
newtype = 'same';
foo = [];

sz = size(M);
dimension = length(sz);

uni = unique(M);
uni(uni==0)=[]; % remove background
n = length(uni);

if max(max(max(M)))>n % Gap in the numerotation
    if dimension == 2
        sz(3)=1;
    end

    % Slow
    if max(max(max(M))) > 1e9
        Mtmp = zeros(sz);
        for k=1:length(uni)
            Mtmp(M==uni(k))=k;
        end
        M=Mtmp;
        clear Mtmp;
    else
        % Faster but can cause memory error if max(max(max(M))) is too high
        conv_uni = zeros(max(uni),1);
        conv_uni(1:uni(1)) = 1;
        for k=2:1:n
            x0 = uni(k-1)+1;
            x1 = uni(k);
            conv_uni(x0:x1) = ones(uni(k)-uni(k-1),1)*k;
        end
        conv_uni = [0; conv_uni]; % Index starts at 1, not 0, in MATLAB

        for z=1:1:sz(3)
            for x=1:1:sz(1)
                for y=1:1:sz(2)
                    M(x,y,z) = conv_uni(M(x,y,z)+1);
                end
            end
        end
    end

end

% Convert to uint8 or uint16
[M] = fct_intconvert(M);

end