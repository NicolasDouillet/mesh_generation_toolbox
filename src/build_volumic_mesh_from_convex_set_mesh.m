function [V_out, T_out] = build_volumic_mesh_from_convex_set_mesh(V_in, T_in)
% build_volumic_mesh_from_convex_set_mesh


C = mean(V_in,1) ; % point set isobarycentre
V_out = cat(1,V_in,C);
C_id = size(V_out,1);

new_tgl = create_vol_mesh_triangles(T_in,C_id);
T_out = cat(1,T_in,new_tgl);


end % build_volumic_mesh_from_convex_set_mesh


% create_vol_mesh_triangles subfunction
function [new_tgl] = create_vol_mesh_triangles(T, C_id)
% create_vol_mesh_triangles : function to create new triangles
% part of tetrahedric volumic mesh (inside triangles)

edg_list = unique(query_edges_list(T,'sorted'),'rows');
new_tgl = cat(2,edg_list,C_id*ones(size(edg_list,1),1));


end % create_vol_mesh_triangles