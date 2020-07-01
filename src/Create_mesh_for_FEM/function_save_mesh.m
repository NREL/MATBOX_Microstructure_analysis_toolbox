function [outputArg1,outputArg2] = function_save_mesh(node_region, elem_region, str_region_name_sav, options)

fullpath=[options.save.mainfolder options.save.stlfolder];
% Export STL
if options.save.microstructurestl==1
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    savestl(node_region,elem_region,[fullpath,'mesh_' str_region_name_sav '_stl.stl'],str_region_name_sav)
end
if options.save.microstructurebinarystl==1
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    savebinstl(node_region,elem_region,[fullpath,'mesh_' str_region_name_sav '_binarystl.stl'],str_region_name_sav)
end

% % Export verteces (nodes) and cells (tetrahedron) - so that you recreate
% manually the mesh in other software easily (such as FEniCS)
fullpath=[options.save.mainfolder options.save.meshdatafolder];
if options.save.microstructuremesh==1
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    elem_=elem_region(:,1:4); % Remove useless last column
    elem_=elem_-1; % Start at 0 for python index convention
    save([fullpath,['Nodes_' str_region_name_sav '.mat']],'node_region')
    save([fullpath,['Tetrahedron_' str_region_name_sav '.mat']],'elem_')
end

end

