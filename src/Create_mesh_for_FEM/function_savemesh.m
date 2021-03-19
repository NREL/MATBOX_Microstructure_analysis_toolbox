function [] = function_savemesh(node_, face_, elem_, subdomain_, str_name, options)

% Export STL
if options.save_stl
    savestl(node_,elem_,[options.folder,'mesh_' str_name '_stl.stl'],str_name)
end
if options.save_binarystl
    savebinstl(node_,elem_,[options.folder,'mesh_' str_name '_binarystl.stl'],str_name)
end

% Export msh
if options.save_msh
    savemsh(node_,elem_(:,1:4),[options.folder,'mesh_' str_name '.msh'])
end

% Export Abaqus
if options.save_inp
    saveabaqus(node_,face_(:,1:3),elem_(:,1:4),[options.folder,'mesh_' str_name '.inp'])
end

% Export verteces (nodes), cells (tetrahedron) and subodmains, so that you can recreate
% manually the mesh in other software easily (such as FEniCS)
elem_=elem_(:,1:4); % Remove useless last column
if options.indexstart_zero
    elem_=elem_-1;
end
if options.save_mat
    save([options.folder,['Nodes_' str_name '.mat']],'node_', '-v7.3')
    save([options.folder,['Tetrahedron_' str_name '.mat']],'elem_', '-v7.3')
    save([options.folder,['Subdomain_' str_name '.mat']],'subdomain_', '-v7.3')
end
if options.save_csv
    writematrix(node_,[options.folder 'Nodes_' str_name '.csv']);
    writematrix(elem_,[options.folder 'Tetrahedron_' str_name '.csv']);
    writematrix(subdomain_,[options.folder 'Subdomain_' str_name '.csv']);
end

% Save min-max coordinates
coor_min_max = [min(node_(:,1)) min(node_(:,2)) min(node_(:,3)); max(node_(:,1)) max(node_(:,2)) max(node_(:,3))];
save([options.folder str_name '_coor_min_max.txt'],'coor_min_max','-ascii');
%dlmwrite([options.folder str_name '_coor_min_max.txt'],coor_min_max);

end
