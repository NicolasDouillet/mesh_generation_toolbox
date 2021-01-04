% test_build_volumic_mesh_from_convex_set_mesh

clear all, close all, clc;

addpath('../src');
addpath('../data');

% Generate the mesh of a convex surface
id = 4;
nb_it = 5;
projection_mode = 'edge_oversamples';
volumic_mesh = false;
[V,T] = build_geoid(id,projection_mode,nb_it,volumic_mesh);
[V,T] = build_volumic_mesh_from_convex_set_mesh(V,T);

% Remove some vertices to have an insight view
V_set = 1:16;
[V,T] = remove_vertices(V_set,V,T);
plot_mesh(V,T);
view(-114,11);