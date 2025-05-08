function [Mseg] = Segm_InstanceRecombineParts(idlargestaxe,bounds_independant,bounds_overlapping,Ind,Over)

nfile = length(Ind) + length(Over);
kind = 0;
kover = 0;

for k = 1:1:nfile
    if rem(k,2)==0
        kover = kover+1;
        inst = Over(kover).Mseg;
    else
        kind = kind+1;
        inst = Ind(kind).Mseg;
    end
    
    if idlargestaxe==1
        p.action = 'Swap axis 1 with axis 3';
        [inst,~,~] = fct_flipswap(inst,p);
    elseif idlargestaxe==2
        p.action = 'Swap axis 2 with axis 3';
        [inst,~,~] = fct_flipswap(inst,p);        
    end
    sz = size(inst);

    if rem(k,2)==1 % Odd number. Remove all labels at the edges
        % Identify id of instances at edges
        if k==1
            outside_ids = unique( inst(:,:,end) );
        elseif k==nfile
            outside_ids = unique( inst(:,:,1) );
        else
            outside_0 = unique( inst(:,:,end) );
            outside_1 = unique( inst(:,:,1) );
            outside_ids = unique([outside_0; outside_1]);
        end

        % Background
        outside_ids(outside_ids==0)=[];

        % Remove such instances
        if ~isempty(outside_ids)
            for kid=1:1:length(outside_ids)
                inst(inst==outside_ids(kid))=0;
            end
        end

    else % Even number. Keep only labels at the middle
        middle_z = round(sz(3)/2);
        inside_ids = unique( inst(:,:,middle_z) );
        BW = zeros(sz);
        for kid = 1:1:length(inside_ids)
            BW(inst==inside_ids(kid))=inside_ids(kid);
        end
        inst = BW;
        [inst] = fct_intconvert(inst); % uint8 or uint16
    end
 
    % Renumerotate
    [inst,~,~] = fct_remvovegapinnumbering(inst,[]);
    [inst] = fct_intconvert(inst);
    if k>1
        id0 = find(inst==0);
        max_id = max(max(max( file(k-1).new_instances )));
        inst = double(inst) + double(max_id);
        inst(id0) = 0;
    end

    file(k).new_instances = inst;
end    

%% RECOMBINE
zmax = bounds_independant(end,end);
Mseg = zeros(sz(1),sz(2),zmax);
overlap = zeros(sz(1),sz(2),zmax);

kind = 0;
kover = 0;
for k = 1:1:nfile
    if rem(k,2)==0
        kover = kover+1;
        z0 = bounds_overlapping(kind,1);
        z1 = bounds_overlapping(kind,2);
    else
        kind = kind+1;
        z0 = bounds_independant(kind,1);
        z1 = bounds_independant(kind,2);
    end

    if k>1 % Check for minor overlapping
        BW1 = Mseg(:,:,z0:z1);
        BW2 = double(file(k).new_instances);
        BW = zeros(size(BW1));
        cond1 = BW1>0;
        cond2 = BW2>0;
        BW(cond1.*cond2 == 1)=1;
        overlap(:,:,z0:z1) = BW;
    end

    Mseg(:,:,z0:z1) = Mseg(:,:,z0:z1) + double(file(k).new_instances);

end

[Mseg,~,~] = fct_remvovegapinnumbering(Mseg,[]);

%% FIX OVERLAPP, IF ANY
idoverlapp = find(overlap==1);

tmp = Mseg;
tmp(idoverlapp)=0;
id0 = find(Mseg==0);
[~,idx] = bwdist(tmp);
tmp(id0) = Mseg(idx(id0));

if ~isempty(idoverlapp)
    [~,idx] = bwdist(~overlap);
    Mseg(idoverlapp) = tmp(idx(idoverlapp));
    [Mseg,~,~] = fct_remvovegapinnumbering(Mseg,[]);
end

%% FIX MISSING PARTICLES, IF ANY

for kind = 1:1:length(Over)
    z0 = bounds_overlapping(kind,1);
    z1 = bounds_overlapping(kind,2);

    initial = double(Over(kover).Mseg);
    if idlargestaxe==1
        p.action = 'Swap axis 1 with axis 3';
        [initial,~,~] = fct_flipswap(initial,p);
    elseif idlargestaxe==2
        p.action = 'Swap axis 2 with axis 3';
        [initial,~,~] = fct_flipswap(initial,p);        
    end

    current = Mseg(:,:,z0:z1);

    cond1 = current==0;
    cond2 = initial~=0;
    id = find(cond1.*cond2);
    if ~isempty(id)
        max_ = max(max(max(Mseg)));
        current(id) = initial(id) + double(max_);
        Mseg(:,:,z0:z1) = current;
    end
   
end
[Mseg,~,~] = fct_remvovegapinnumbering(Mseg,[]);

%% RE-ORIENT
if idlargestaxe==1
    p.action = 'Swap axis 1 with axis 3';
    [Mseg,~,~] = fct_flipswap(Mseg,p);
elseif idlargestaxe==2
    p.action = 'Swap axis 2 with axis 3';
    [Mseg,~,~] = fct_flipswap(Mseg,p);
end

end