function [facets] = Function_calculate_facets(elem)

% given a reference tet with edges 1,2,3,4, list the faces
facet_pattern = nchoosek([1,2,3,4],3);

% list all facets for every test in elem
facets = elem(:,facet_pattern);

% process to get the right shape
facets = reshape(facets,[4*length(elem),3]);

% sort each row so unique will work
facets = sort(facets,2);

% eleminate extra facets
facets = unique(facets,'rows');

% sort the rows to match fenics
facets = sortrows(facets,1);

end

