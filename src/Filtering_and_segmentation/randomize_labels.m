function [Label_rnd,n] = randomize_labels(Label)

checkdist = true;
d_check = 100;

uni = unique(Label);
uni(uni==0)=[]; % remove background
n = length(uni);

if ~checkdist
    rnd = randperm(n); % random permutation of the integers from 1 to n without repeating elements.
else
    choices = 1:1:n;
end   

Label_rnd = zeros(size(Label));
for k=1:1:n
    if checkdist
        sz = size(Label);
        dimension = length(sz);
        idx = find(Label==uni(k));
        if dimension == 2
            [IX, IY] = ind2sub(sz,idx);
            x_min = max([1 min(IX)-d_check]); x_max = min([sz(1) max(IX)+d_check]);
            y_min = max([1 min(IY)-d_check]); y_max = min([sz(2) max(IY)+d_check]);
            sub = Label_rnd(x_min:x_max, y_min:y_max);
        else
            [IX, IY, IZ] = ind2sub(sz,idx);
            x_min = max([1 min(IX)-d_check]); x_max = min([sz(1) max(IX)+d_check]);
            y_min = max([1 min(IY)-d_check]); y_max = min([sz(2) max(IY)+d_check]);
            z_min = max([1 min(IZ)-d_check]); z_max = min([sz(3) max(IZ)+d_check]);
            sub = Label_rnd(x_min:x_max, y_min:y_max, z_min:z_max);
        end
        unisub = unique(sub);
        unisub(unisub==0)=[];
        n_sunisub = length(unisub);
        n_choice = length(choices);
        if n_sunisub==0
            r = randi(n_choice);
            Label_rnd(idx) = choices(r);
            choices(r) = [];
        else
            d = zeros(n_sunisub,n_choice);
            for kk=1:1:n_sunisub
                d(kk,:) = abs(choices - unisub(kk));
            end
            m = min(d);
            id = find(m==max(m));
            Label_rnd(idx) = choices(id(1));
            choices(id(1)) = [];
        end
    else
        Label_rnd(Label==uni(k)) = rnd(k);
    end
end
end




