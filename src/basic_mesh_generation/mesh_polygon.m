function [V, T] = mesh_polygon(P, mode)
%% mesh_polygon : function to triangularly mesh a given polygon.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2021-2025.
%
%
%%% Input arguments
%
%       [ |  |  |]
% - P = [Px Py Pz], real matrix double, the polygon, size(P) = [nb_vertices,3]. Mandatory.
%       [ |  | | ]
%
% - mode : character string in the set {'raw','sorted'},
%          the type of contour considered. Case insensitive. Mandatory.
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
% P at least must be a quadrangle; Function won't perform on polygons which already are triangles.
% Triangles / normals may not be coherently oriented.


%% Body
[V,T] = discrete_contour_mesh_patch(P,mode);


end % mesh_polygon