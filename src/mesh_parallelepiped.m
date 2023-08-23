function [V, T] = mesh_parallelepiped(L, d, h, edg_nb_smpl)
%
% Author & support : nicolas.douillet (at) free.fr, 2023.


% L : length (X)
% d : depth  (Y)
% h : height (Z)

% Parallelepiped is centred on the origin, [0 0 0].
%
% Normals not uniformally oriented here.


[V,T] = mesh_cube(edg_nb_smpl);

V(:,1) = L*V(:,1);
V(:,2) = d*V(:,2);
V(:,3) = h*V(:,3);


end % mesh_rectangle_parallelepiped