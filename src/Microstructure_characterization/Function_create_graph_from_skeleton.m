function [G] = Function_create_graph_from_skeleton(node_coordinate,node_connectivity)

%% INPUT DATA

% Variable: node_coordinate
% Particle_id (integer)
% mass center coordinate X (micrometer)
% mass center coordinate Y (micrometer)
% mass center coordinate Z (micrometer)
% Volume of the particle (number of voxels)
% Equivalent diameter (micrometer)
% Note: node_coordinate(i,1)=i

% Variable: node_connectivity
% node_connectivity(i,j)=number of face that connect two particles
% The particle id of these two particles is node_connectivity(i,1) and node_connectivity(1,j)
% node_connectivity is symmetric
% node_connectivity(0,0)=1, node_connectivity(i,i)=1 for i>1


%% GRAPH SYNTAX

% The syntax is (for undirected graph):
% node_A     = [1 1 1 2 2 3 3 4 5 5 6 7];
% node_B     = [2 4 8 3 7 4 6 5 6 8 7 8];
% weights_AB = [10 10 1 10 1 10 1 1 12 12 12 12];
% G = graph(node_A,node_B,weights_AB)


%% GET NODES

% Initialise nodes and weights
node_A=[];
node_B=[];
weights_AB=[];

node_A_bis=[];
node_B_bis=[];
weights_AB_bis=[];

% Number of particles;
number_particle = length(node_coordinate(:,1));

% Create node to node connection by using the node_connectivity matrix
for current_particle=1:1:number_particle
    % Get node A
    current_node_A=node_connectivity(current_particle+1,1);
    % Get all node B connected to node A
    index_ = find(node_connectivity(current_particle+1,:)~=0);
    
    % % The following lines lead to duplicate all connections
    % Substracte 1 
    % index_=index_-1;
    % Remove 0 and current_node_A
    % index_(index_==0)=[];
    % index_(index_==current_node_A)=[];
    
    % Instead of, we will consider only the upper-right of the symmetric node_connectivity matrix
    % Substracte 1 
    index_=index_-1;
    % Kepp only the upper-right
    index_(index_<=current_node_A)=[];
    % Number of connection
    connection_number = length(index_);    
    if connection_number>0
        
        % % Update node-node connection for the graph
        % Update node_A
        node_A=[node_A ones(1,connection_number)*current_node_A];
        
        % Update node_B
        node_B=[node_B index_];
        
        % % Update node-node weight connection for the graph
        % Get coordinate of node A
        x_A = ones(connection_number,1)*node_coordinate(current_particle,2);
        y_A = ones(connection_number,1)*node_coordinate(current_particle,3);
        z_A = ones(connection_number,1)*node_coordinate(current_particle,4);
        % Get ALL coordinate of nodes B
        x_B = node_coordinate(index_,2);
        y_B = node_coordinate(index_,3);
        z_B = node_coordinate(index_,4);
        % Calculate euclidean distances
        euclidean_distance = ((x_B-x_A).^2 + (y_B-y_A).^2 +(z_B-z_A).^2).^(0.5);
        % Update weights_AB
        weights_AB=[weights_AB euclidean_distance'];
        
    else
        node_A=[node_A current_node_A];
        node_B=[node_B current_node_A];    
        weights_AB=[weights_AB 0];
end

%% CREATE THE GRAPH
G = graph(node_A,node_B,weights_AB);

%plot(G);


end