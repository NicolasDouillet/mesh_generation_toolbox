function [V, T] = mesh_quadrangle(M1, M2, M3, M4, edg_nb_smpl)
%
% Author & support : nicolas.douillet (at) free.fr, 2023.


[V1, T1] = sample_triangle(M1',M2',M3',edg_nb_smpl);
[V2, T2] = sample_triangle(M3',M4',M1',edg_nb_smpl);

V = cat(1,V1,V2);
T = cat(1,T1,T2+size(V1,1));


end % mesh_quadrangle