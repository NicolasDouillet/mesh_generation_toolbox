function [V_out, T_out] = build_volumic_mesh_from_convex_set_mesh(V_in, T_in)
%% build_volumic_mesh_from_convex_set_mesh : function to build a volumic mesh from a given input convex surface mesh.
% 
%
% Author & support : nicolas.douillet (at) free.fr, 2021.
%
%
% Input arguments
%
%          [ |    |    |  ]
% - V_in = [X_in Y_in Z_in], real matrix double, the input point set, size(V_in) = [nb_input_vertices,3].
%          [ |    |    |  ]
%
%          [  |     |     |  ]
% - T_in = [i1_in i2_in i3_in], positive integer matrix double, the input triangulation, size(T_in) = [nb_input_triangles,3].
%          [  |     |     |  ]
%
%
% Output arguments
%
%           [  |     |     |  ]
% - V_out = [X_out Y_out Z_out], real matrix double, the output point set, size(V_out) = [nb_output_vertices,3],
%           [  |     |     |  ]
%
%           where nb_output_vertices = nb_input_vertices - nb_duplicata.
%
%           [  |      |      |   ]
% - T_out = [i1_out i2_out i3_out], positive integer matrix double, the output triangulation, size(T_out) = [nb_output_triangles,3].
%           [  |      |      |   ]

%% Body
C = mean(V_in,1); % point set isobarycentre
V_out = cat(1,V_in,C);
C_id = size(V_out,1);

new_tgl = create_vol_mesh_triangles(T_in,C_id);
T_out = cat(1,T_in,new_tgl);


end % build_volumic_mesh_from_convex_set_mesh


%% create_vol_mesh_triangles subfunction
function [new_tgl] = create_vol_mesh_triangles(T, C_id)
% create_vol_mesh_triangles : function to create new triangles
% part of tetrahedric volumic mesh (inside triangles)

edg_list = unique(query_edges_list(T,'sorted'),'rows');
new_tgl = cat(2,edg_list,C_id*ones(size(edg_list,1),1));


end % create_vol_mesh_triangles