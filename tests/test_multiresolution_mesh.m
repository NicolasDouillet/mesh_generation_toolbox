% test multiresolution_mesh

clear all, close all, clc;

addpath('../src');
addpath('../data');


% nb_vtx = 128;
% X = 2*(rand(nb_vtx,1)-0.5);
% Y = 2*(rand(nb_vtx,1)-0.5);
% Z = 2*(rand(nb_vtx,1)-0.5);
% 
% Rho = X.^2 + Y.^2 + Z.^2;
% i = Rho <= 1;
% X = X(i);
% Y = Y(i);
% Z = Z(i);
% V = cat(2,X,Y,Z);


load('concave_Reuleaux_tetrahedron.mat');

T = multiresolution_mesh(V,1);

% plot_point_set(V);
% 


plot_mesh(V,T);
alpha(0.5)


% check self intersections