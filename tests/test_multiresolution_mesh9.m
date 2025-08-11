% test multiresolution_mesh9

% function T_new =  test_multiresolutiomesh9(V, k)

clc;

addpath(genpath('../src/'));
addpath('../data/');

addpath('C:\Users\Nicolas\Desktop\TMW_contributions\mesh_processing_toolbox\data');

close all;

% nb_vtx = 128;
% V = 2*(rand(nb_vtx,3)-0.5);
% V = V./vecnorm(V')';
% V = V*0.1.*(1+rand(nb_vtx,1));


% Tests avec :
%
% - sinusoidal icosahedron / dodecahedron
% - visage
% - logo Matlab

% load('spiky_cell_like_surface.mat'); % beware of potential duplicated vertices
% -> à vérifier

% load('sinusoidal_icosahedron_ULR.mat');
load('concave_Reuleaux_tetrahedron.mat');
% load('Archi_spiral.mat');
% load('Gargoyle_5k.mat');
% load('meshed_mtlb_logo.mat');
clear T;

% V = unique(V,'rows'); % if necessary (presence of duplicated vertices)

% k = 3;
% [V,T] = multiresolution_mesh9(V,k);
[V,T] = multiresolution_mesh11(V);
select_face_normals(V,T);
% plot_point_set_and_mesh(V,T);
% plot_point_set(V);

% N = fliplr(face_normals(V,T));
% E = query_edg_list(T,'sorted');
% select_edge_normals(V,T,N,E);


alpha(1);
camlight left;