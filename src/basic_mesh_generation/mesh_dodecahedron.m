function [V, T] = mesh_dodecahedron(edg_nb_smpl)
%% mesh_dodecahedron : function to mesh a dodecahedron.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2023-2025.
%
%
%%% Input argument
%
% - edg_nb_smpl : positive integer scalar double, the number of samples for each one the dodecahedron edges. Mandatory.
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
% Dodecahedron is centered on the origin, [0 0 0].
% Triangles / normals are coherently oriented and facing outward.


%% Body
[M,F] = platonic_solids(5,1,'triangle'); % dodecahedron made of 12 triangles at this point

V = zeros(0,3);
T = zeros(0,3);


for k = 1:size(F,1) % 12
    
    [Vk,Tk] = mesh_triangle(M(F(k,3),:)',M(F(k,2),:)',M(F(k,1),:)',edg_nb_smpl);        
        
    T = cat(1,T,Tk+size(V,1));
    V = cat(1,V,Vk);

end


end % mesh_dodecahedron