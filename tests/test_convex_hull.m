% test convex_hull

clear all, close all, clc;

addpath('../src');
addpath('../data');


nb_vtx = 128;
X = 2*(rand(nb_vtx,1)-0.5);
Y = 2*(rand(nb_vtx,1)-0.5);
Z = 2*(rand(nb_vtx,1)-0.5);

Rho = X.^2 + Y.^2 + Z.^2;
i = Rho <= 1;
X = X(i);
Y = Y(i);
Z = Z(i);
V = cat(2,X,Y,Z);

% filenames = {'concave_Reuleaux_tetrahedron';...
%             };


% id = 1;
% filename = strcat(cell2mat(filenames(id,1)),'.mat');         
% load(filename);


T_cv = convhull(V(:,1),V(:,2),V(:,3)); % Matlab (R) embeded function
[V_out,T_qcv] = quick_hull(V);         % my quickhull algorithm


plot_mesh(V,T_cv);
axis equal;

plot_mesh(V_out,T_qcv);
axis equal;