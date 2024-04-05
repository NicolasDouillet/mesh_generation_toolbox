% test build_triangulation_from_edge_list


clear all, close all, clc;

addpath('../src');
addpath('../data/');


filenames = {'concave_Reuleaux_tetrahedron'};

filename = strcat(cell2mat(filenames(1,1)),'.mat');         
load(filename);

% Randomly mess up the edge set
E = query_edges_list(T);
new_idx = randperm(size(E,1));
E = E(new_idx,:);
E = fliplr(E);

clear T;

% Build the mesh
T = build_triangulation_from_edge_list(E);
plot_mesh(V,T);