function T_out = switch_neighbor_pair_triangles(T_in, tid1, tid2)
%% switch_neighbor_pair_triangles : function to swuitch
% neighbor pair triangles.
%
% Author : nicolas.douillet9 (at) gmail.com, 2024.
%
%
% Input arguments :
%
%          [  |     |     |  ]
% - T_in = [i1_in i2_in i3_in], positive integer matrix double, the input triangulation, size(T_in) = [nb_input_triangles,3].
%          [  |     |     |  ]
%
% - tid1 : integer scalar double, the index of the first triangle of the pair to switch. 1 <= tid1 <= size(T_in,1).
%
% - tid2 : integer scalar double, the index of the second triangle of the pair to switch. 1 <= tid2 <= size(T_in,1).
%
%
% Output argument :
%
%           [  |      |      |   ]
% - T_out = [i1_out i2_out i3_out], positive integer matrix double, the output triangulation, size(T_out) = [nb_output_triangles,3].
%           [  |      |      |   ]
%
%
% About / others information
%
% May mess up triangle normals orientation


%% Body
% Find common edge
cmn_edge = intersect(T_in(tid1,:),T_in(tid2,:));

% Create new triangles
new_tgl1 = cat(2,setdiff(T_in(tid1,:),cmn_edge),cmn_edge(1),setdiff(T_in(tid2,:),cmn_edge));
new_tgl2 = cat(2,setdiff(T_in(tid2,:),cmn_edge),cmn_edge(2),setdiff(T_in(tid1,:),cmn_edge));

% Update triangulation
T_out = add_triangles(new_tgl1,T_in);
T_out = add_triangles(new_tgl2,T_out);
T_out = remove_triangles(T_out,[tid1,tid2]);


end % switch_neighbor_pair_triangles