function [V, T] = mesh_cube(edg_nb_smpl)
%
% Author & support : nicolas.douillet (at) free.fr, 2023.

% edg_nb_smpl is the number of samples for each one the cube edge.
%
% Cube is centred on the origin, [0 0 0].
%
% Normals not uniformally oriented here.





[M,F] = mesh_platonic_solids(2,1,false,'triangle'); % cube made of 12 triangles at this point

V = zeros(0,3);
T = zeros(0,3);


for k = 1:size(F,1) % 12
    
    [Vk,Tk] = sample_triangle(M(F(k,1),:)',M(F(k,2),:)',M(F(k,3),:)',edg_nb_smpl);        
        
    T = cat(1,T,Tk+size(V,1));
    V = cat(1,V,Vk);

end


end % mesh_cube