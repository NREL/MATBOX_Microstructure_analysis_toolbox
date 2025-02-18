f = ones(50,50,100);

f(1:5,:,:) = 0;
f(10:20,10:20,20:50) = 2;
f(30:40,10:20,40:70) = 3;
f(10:20,10:20,60:100) = 4;
file(1).instances = f;
file(1).z0 = 1;
file(1).z1 = 100;
file(2).instances = f;
file(2).z0 = 51;
file(2).z1 = 150;
file(3).instances = f;
file(3).z0 = 101;
file(3).z1 = 200;

% clearvars
% close all
% clc

%% IMPORT FILES
has_background = true;
id_background = 1;


%% PER FILE
nfile = length(file);
for k = 1:1:nfile
    inst = file(k).instances;
    sz = size(inst);

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
    if has_background
        outside_ids(outside_ids==id_background)=[];   
    end

    % Remove such instances
    if ~isempty(outside_ids)
        for kid=1:1:length(outside_ids)
            if has_background
                inst(inst==outside_ids(kid))=id_background;
            else
                inst(inst==outside_ids(kid))=0;
            end
        end
    end    

    % Renumerotate
    [inst,~,~] = fct_remvovegapinnumbering(inst,[]);
    if k>1
        id0 = find(inst==0);
        if has_background
            idb = find(inst==id_background);
        end
        max_id = max(max(max( file(k-1).new_instances )));
        inst = inst + max_id;
        inst(id0) = 0;
        if has_background
            inst(idb) = id_background;
        end
    end

    file(k).new_instances = inst;
    %file(k).new_instances = inst(:,:,1:end-5);
end

%% RECOMBINE
res = zeros(sz(1),sz(2),file(nfile).z1) - 1 ;
for k = 1:1:nfile
    sz = size(file(k).new_instances);
    z1 = file(k).z0 + sz(3) - 1;
    res(:,:,file(k).z0:z1) = file(k).new_instances;
end
[res,~,~] = fct_remvovegapinnumbering(res,[]);

%% FILL GAP IN ANY
idunassigned = find(res==-1);
if ~isempty(idunassigned)
    BW = zeros(size(res));
    BW(idunassigned) = 1;
    [~,idx] = bwdist(~BW);
    res(idunassigned) = res(idx(idunassigned));
end

%Microstructure_basic_visualization_interface(res)