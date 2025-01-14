function [V, T] = mesh_quadrangle(M1, M2, M3, M4, edg_nb_smpl)
%% mesh_quadrangle : function to mesh the (M1M2M3M4) quadrangle.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2023-2025.
%
%
%%% Input argument
%
% - edg_nb_smpl : positive integer scalar double, the number of samples. Mandatory.
%
%
%%% Output arguments
%
%        [| | |]
% - V_ = [X Y Z], real matrix double, the output point set, size(V) = [nb_vertices,3]
%        [| | |]
%                 with nb_vertices is the same number as the number of vertices in P.
%
%       [|  |  |]
% - T = [i1 i2 i3], positive integer matrix double, the triangulation, size(T) = [nb_triangles,3].
%       [|  |  |]
%
%
%%% About / other information
%
% Triangles / normals are coherently oriented.


%% Body
[V1, T1] = mesh_triangle(M1',M2',M3',edg_nb_smpl);
[V2, T2] = mesh_triangle(M3',M4',M1',edg_nb_smpl);

V = cat(1,V1,V2);
T = cat(1,T1,T2+size(V1,1));


end % mesh_quadrangle