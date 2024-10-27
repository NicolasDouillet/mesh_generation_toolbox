function [V, T] = volumic_mesh_tetrahedron(nb_it, V1, V2, V3, V4)
%% volumic_mesh_tetrahedron : function to build a volumic mesh of a given
% tetrahedron. Performs barycentre based subdivisions.
%
% Author : nicolas.douillet (at) free.fr, 2024.
%
%
% Input argument
%
% - nb_it : positive integer scalar double, the number of iteration.
%           Optional. Set to 1 by default.
%
% - V1 = [V1x V1y V1z], real row vector double, the first  vertex of the tetrahedron, size(V1) = [1,3]. Optional.
%
% - V2 = [V2x V2y V2z], real row vector double, the second vertex of the tetrahedron, size(V2) = [1,3]. Optional.
%
% - V3 = [V3x V3y V3z], real row vector double, the third  vertex of the tetrahedron, size(V3) = [1,3]. Optional.
%
% - V4 = [V4x V4y V4z], real row vector double, the fourth vertex of the tetrahedron, size(V4) = [1,3]. Optional.
%
%
% Output arguments
%
%       [|  |  | ]
% - V = [Vx Vy Vz], real matrix double, the output point set, size(V) = [nb_output_vertices,3],
%       [|  |  | ]
%
%       [|  |  | ]
% - T = [t1 t2 t3], positive integer matrix double, the triangulation, size(T) = [nb_triangles,3].
%       [|  |  | ]


%% Input parsing
if nargin < 5
    
    % Default basic tetrahedron included in the unit sphere
    V1 = [2*sqrt(2)/3 0 -1/3];
    V2 = [-sqrt(2)/3 sqrt(6)/3 -1/3];
    V3 = [-sqrt(2)/3 -sqrt(6)/3 -1/3];
    V4 = [0 0 1];
    
    if nargin < 1
        
        nb_it = 1;
        
    end
    
end


%% Body
V = cat(1,V1,V2,V3,V4);

T_124 = [1 2 4];
T_234 = [2 3 4];
T_314 = [3 1 4];
T_132 = [1 3 2];

T = cat(1,T_124,T_234,T_314,T_132);
tetra_list = [1 2 3 4];
tetra_start_idx = 1;

for j = 1:nb_it
    
    nb_tetra = size(tetra_list,1);
    
    for k = tetra_start_idx:nb_tetra
        
        [V,T,new_tetra] = split_tetrahedron(V,T,tetra_list(k,:));
        tetra_list = cat(1,tetra_list,new_tetra);
        
    end
    
    tetra_start_idx = nb_tetra + 1;
    
end

T = remove_duplicated_triangles(T);


end % volumic_mesh_tetrahedron


%% split_tetrahedron subfunction
function [V, T, new_tetra] = split_tetrahedron(V, T, tetra_list_k) 
%
% Author : nicolas.douillet (at) free.fr, 2024.


new_vtx = mean(V(tetra_list_k,:),1);
V = cat(1,V,new_vtx);
new_vtx_idx = size(V,1);

new_tgl = cat(2,combnk(tetra_list_k,2),new_vtx_idx*ones(6,1));
T = cat(1,T,new_tgl);

new_tetra = cat(2,combnk(tetra_list_k,3),new_vtx_idx*ones(4,1));


end % split_tetrahedron