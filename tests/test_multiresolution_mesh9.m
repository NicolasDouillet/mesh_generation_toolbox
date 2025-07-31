% test multiresolution_mesh9

% function T_new =  test_multiresolutiomesh9(V, k)

clc;

addpath(genpath('../src/'));
addpath('../data/');


% nb_vtx = 128;
% V = 2*(rand(nb_vtx,3)-0.5);
% V = V./vecnorm(V')';
% V = V*0.1.*(1+rand(nb_vtx,1));


% Tests avec :
%
% - surface boite à oeufs
% - surface spiralante,
% - visage
% - logo Matlab

load('concave_Reuleaux_tetrahedron.mat');
% load('Gargoyle_5k.mat');
% load('meshed_mtlb_logo.mat');

% V = unique(V,'rows'); % if necessary (presence of duplicated vertices)

k = 6;
[~,T_new] = multiresolution_mesh9(V,k);
plot_mesh(V,T_new);
alpha(1);
camlight left;