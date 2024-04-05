function [V, T] = mesh_octahedron(edg_nb_smpl)
% mesh_octahedron : function to mesh a octahedron.
%
% Author : nicolas.douillet (at) free.fr, 2023-2024.
%
%
%%% Input arguments :
%
% - edg_nb_smpl : positive integer scalar double, the number of samples for each one the octahedron edges.
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
%%% About / others information
%
% Octahedron is centered on the origin, [0 0 0].
% Triangles / normals are coherently oriented and facing outward.


[M,F] = platonic_solids(3,1); % octahedron made of 8 triangles at this point

V = zeros(0,3);
T = zeros(0,3);


for k = 1:size(F,1) % 8
    
    [Vk,Tk] = mesh_triangle(M(F(k,3),:)',M(F(k,2),:)',M(F(k,1),:)',edg_nb_smpl);        
        
    T = cat(1,T,Tk+size(V,1));
    V = cat(1,V,Vk);

end


end % mesh_octahedron