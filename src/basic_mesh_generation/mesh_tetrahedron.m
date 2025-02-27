function [V, T] = mesh_tetrahedron(edg_nb_smpl)
%% mesh_tetrahedron : function to mesh a tetrahedron.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2023-2025.
%
%
%%% Input argument
%
% - edg_nb_smpl : positive integer scalar double, the number of samples for each one the tetrahedron edges. Mandatory.
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
%%% About / others information
%
% Tetrahedron is centered on the origin, [0 0 0].
% Triangles / normals are coherently oriented and facing outward.


%% Body
[M,F] = platonic_solids(1,1); % tetrahedron made of 4 triangles at this point

V = zeros(0,3);
T = zeros(0,3);


for k = 1:size(F,1) % 4
    
    [Vk,Tk] = mesh_triangle(M(F(k,3),:)',M(F(k,2),:)',M(F(k,1),:)',edg_nb_smpl);        
        
    T = cat(1,T,Tk+size(V,1));
    V = cat(1,V,Vk);

end


end % mesh_tetrahedron