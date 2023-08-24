function [V, T] = mesh_ellipse(a, b, edg_nb_smpl)
% mesh_ellipse : function to mesh an ellipse, based on the disk.
%
% Author & support : nicolas.douillet (at) free.fr, 2023.
%
%
%%% Input arguments :
%
% - a : positive real scalar double, X axis multiplicating factor.
%
% - b : positive real scalar double, Y axis multiplicating factor.
%
% - edg_nb_smpl : positive integer scalar double, one sixth of the number of samples on the ellipse perimeter.
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
% Ellipse is centered on the origin, [0 0 0], and has (Oz) for normal axis.
% Normals are coherently oriented.


[V,T] = mesh_disk(1,edg_nb_smpl);

V(:,1) = a * V(:,1);
V(:,2) = b * V(:,2);


end % mesh_ellipsoid