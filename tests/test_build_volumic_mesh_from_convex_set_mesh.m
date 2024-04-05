% test_volumic_mesh_from_convex_set_mesh


clear all, close all, clc;

addpath('../src');
addpath('../data');

% Generate the mesh of a convex surface
id = 4;
nb_it = 5;
sampling_mode = 'edge';
volumic_mesh = false;
[V,T] = mesh_geoid(id,nb_it,sampling_mode);
[V,T] = volumic_mesh_from_convex_set_mesh(V,T);

plot_mesh(V,T);