% test dual_of_trimesh


clear all, close all, clc;

addpath('../src');
addpath('../data');


filenames = {'concave_Reuleaux_tetrahedron'};

         
filename = strcat(cell2mat(filenames(1,1)),'.mat');         
load(filename);

plot_mesh(V,T);

[V_dual,T_dual] = dual_of_trimesh(V,T);
plot_mesh(V_dual,T_dual);       