function [V, T] = mesh_parallelepiped(L, d, h, edg_nb_smpl)
% mesh_parallelepiped : function to mesh a parallelepiped.
%
% Author : nicotangente (at) free.fr, 2023-2024.
%
%
%%% Input arguments :
%
% - L : real scalar double, the length (X) of the parallelepiped.
% - d : real scalar double, the depth  (Y) of the parallelepiped.
% - h : real scalar double, the height (Z) of the parallelepiped.
%
%
%%% Output arguments :
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
%%% About / other informations
%
% Parallelepiped is centered on the origin, [0 0 0].
% Triangles / normals are coherently oriented and facing outward.


[V,T] = mesh_cube(edg_nb_smpl);

V(:,1) = L*V(:,1);
V(:,2) = d*V(:,2);
V(:,3) = h*V(:,3);


end % mesh_rectangle_parallelepiped