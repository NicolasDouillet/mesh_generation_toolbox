function [V, T] = mesh_cube(edg_nb_smpl)
%% mesh_cube : function to mesh a cube.
%
% Author : nicolas.douillet9 (at) gmail.com, 2023-2024.
%
%
% Input arguments :
%
% - edg_nb_smpl : positive integer scalar double, the number of samples for each one the cube edges.
%
%
% Output arguments :
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
% About / others information
%
% Cube is centered on the origin, [0 0 0].
% Triangles / normals are coherently oriented and facing outward.


%% Body
[M,F] = platonic_solids(2,1,'triangle'); % cube made of 12 triangles at this point

V = zeros(0,3);
T = zeros(0,3);


for k = 1:size(F,1) % 12
    
    [Vk,Tk] = mesh_triangle(M(F(k,3),:)',M(F(k,2),:)',M(F(k,1),:)',edg_nb_smpl);        
        
    T = cat(1,T,Tk+size(V,1));
    V = cat(1,V,Vk);

end


end % mesh_cube