% test build_triangulation_from_edge_list

clc;

addpath(genpath('../src'));
addpath('../data');


filenames = {'concave_Reuleaux_tetrahedron'};

filename = strcat(cell2mat(filenames(1,1)),'.mat');         
load(filename);

% Randomly mess up the edge set
E = query_edges_list(T);
new_id = randperm(size(E,1));
E = E(new_id,:);
E = fliplr(E);

clear T;

% Build the mesh
T = build_triangulation_from_edge_list(E);
select_face_normals(V,T);