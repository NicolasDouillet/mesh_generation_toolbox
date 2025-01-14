function [V, T] = mesh_ellipsoid(a, b, c, nb_samples)
%% mesh_ellipsoid : function to mesh an ellipsoid based on the geoid,
% based itself on a subsampling of the icosahedron.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2023-2025.
%
%
%%% Input arguments
%
% - a : positive real scalar double, X axis multiplicating factor. Mandatory.
%
% - b : positive real scalar double, Y axis multiplicating factor. Mandatory.
%
% - c : positive real scalar double, Z axis multiplicating factor. Mandatory.
%
% - nb_samples : positive integer scalar, the number of samples used to
%                subsample each triangle of the icosahedron. Mandatory.
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
% Ellipsoid is centered on the origin, [0 0 0].
% Triangles / normals are coherently oriented and facing outward.


%% Body
[V,T] = mesh_geoid(4,nb_samples,'edge');

V(:,1) = a * V(:,1);
V(:,2) = b * V(:,2);
V(:,3) = c * V(:,3);


end % mesh_ellipsoid