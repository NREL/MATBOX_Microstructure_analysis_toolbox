function [M_rnd,newtype,foo] = fct_randomize(M,p)
newtype = 'same';
foo = [];

sz = size(M);
dimension = length(sz);

checkdist = p.basedondistance;
d_check = p.dist;

uni = unique(M);

if p.keepunchangedsomelabels
    cst  = str2num(p.cst);
    n_cst = length(cst);
    for k=1:1:n_cst
        uni(uni==cst(k))=[]; % remove constant labels
    end
end
n = length(uni);

M_rnd = zeros(size(M));

% Gap in the numerotation, let's fix this first if needed
[M,~,~] = fct_remvovegapinnumbering(M,[]);

if checkdist

    % choices = 1:1:n;
    % for k=1:1:n
    %     idx = find(M==uni(k)); %% THIS IS TOO SLOW
    %     if dimension == 2
    %         [IX, IY] = ind2sub(sz,idx);
    %         x_min = max([1 min(IX)-d_check]); x_max = min([sz(1) max(IX)+d_check]);
    %         y_min = max([1 min(IY)-d_check]); y_max = min([sz(2) max(IY)+d_check]);
    %         sub = M_rnd(x_min:x_max, y_min:y_max);
    %     else
    %         [IX, IY, IZ] = ind2sub(sz,idx);
    %         x_min = max([1 min(IX)-d_check]); x_max = min([sz(1) max(IX)+d_check]);
    %         y_min = max([1 min(IY)-d_check]); y_max = min([sz(2) max(IY)+d_check]);
    %         z_min = max([1 min(IZ)-d_check]); z_max = min([sz(3) max(IZ)+d_check]);
    %         sub = M_rnd(x_min:x_max, y_min:y_max, z_min:z_max);
    %     end
    %     unisub = unique(sub);
    %     unisub(unisub==0)=[];
    %     n_sunisub = length(unisub);
    %     n_choice = length(choices);
    %     if n_sunisub==0
    %         r = randi(n_choice);
    %         M_rnd(idx) = choices(r);
    %         choices(r) = [];
    %     else
    %         d = zeros(n_sunisub,n_choice);
    %         for kk=1:1:n_sunisub
    %             d(kk,:) = abs(choices - unisub(kk));
    %         end
    %         m = min(d);
    %         id = find(m==max(m));
    %         M_rnd(idx) = choices(id(1));
    %         choices(id(1)) = [];
    %     end
    % end   

    choices = 1:1:n;
    %rem_choices = choices;
    seen = [];
    rnd = zeros(n,1);

    idx = find(M~=0);
    if dimension == 2
        [IX,IY] = ind2sub(sz,idx);
    else
        [IX,IY,IZ] = ind2sub(sz,idx);
    end

    for k=1:1:length(idx)
        %k/length(idx)
        val = M(idx(k));
        if isempty(seen)
            rnd(val) = choices(1);
            M_rnd(idx(k)) = choices(1);
            seen = val;
            %rem_choices(1)=[];
        else
            if sum(sum(seen==val)) % Already assigned
                M_rnd(idx(k)) = rnd(val);
            else
                x_min = max([1 IX(k)-d_check]); x_max = min([sz(1) IX(k)+d_check]);
                y_min = max([1 IY(k)-d_check]); y_max = min([sz(2) IY(k)+d_check]);
                if dimension == 2
                    sub=M_rnd(x_min:x_max,y_min:y_max);
                else
                    z_min = max([1 IZ(k)-d_check]); z_max = min([sz(3) IZ(k)+d_check]);
                    sub=M_rnd(x_min:x_max,y_min:y_max,z_min:z_max);
                end
                unisub = unique(sub);
                unisub(unisub==0)=[];

                %if isempty(unisub) % No constraint on choice
                %    rnd(val) = rem_choices(1);
                %    M_rnd(idx(k)) = rem_choices(1);
                %    rem_choices(1)=[];
                %else % Not alone, we need to choose a label far from all others already present
                    nsub = length(unisub);
                    d = zeros(n,nsub);
                    for kk=1:1:nsub
                        d(:,kk) = abs(choices-unisub(kk)); % Calculate distance
                    end
                    d = min(d,[],2); % Minimum distance considering all potential choices
                    idm = find(d==max(d)); % Take maximum distance
                    rnd(val) = choices(idm(1));
                    M_rnd(idx(k)) = rnd(val);
                    %rem_choices(rem_choices==choices(idm(1)))=[];
                    
                %end
                seen = [seen val];

            end
        end
    end


else

    rnd = randperm(n);

    % One liner but so slow...
    %for k=1:1:n
    %    M_rnd(M==uni(k)) = rnd(k);
    %end

    % Dumber (3 for loops!) but so faster
    % rnd = [0 rnd];
    % if dimension == 2
    %     sz(3)=1;
    % end
    % for z=1:1:sz(3)
    %     for x=1:1:sz(1)
    %         for y=1:1:sz(2)
    %             M_rnd(x,y,z) = rnd(M(x,y,z)+1);
    %         end
    %     end
    % end


    if p.keepunchangedsomelabels
        rnd = rnd + max(cst);
        rnd =[cst rnd];
    end
    if dimension == 2
        sz(3)=1;
    end
    for z=1:1:sz(3)
        for x=1:1:sz(1)
            for y=1:1:sz(2)
                M_rnd(x,y,z) = rnd(M(x,y,z)+1);
            end
        end
    end    

end

% Convert to uint8 or uint16
[M_rnd] = fct_intconvert(M_rnd);

end