% test select_edge_normals

clc;

addpath(genpath('../src/'));
addpath('../data');


[V,T] = platonic_solids(5,1,'default');


N = face_normals(V,T);
E = query_edg_list(T,'sorted');
N_e = select_edge_normals(V,T,N,E);