% test switch_neighbor_pair_triangles

clc;

addpath(genpath('../src/'));
addpath('../data/');

I = [1 0 0];
J = [0 1 0];
K = [0 0 1];
O = [0 0 0];

V = cat(1,I,J,K,O);

T = [4 1 2;...
     4 2 3];
 
select_face_normals(V,T); 
T = switch_neighbor_pair_triangles(T,1,2);
select_face_normals(V,T); 