function tgl_id_list = find_triangle_indices_from_edg_list(T, edg_list)
%% find_triangle_indices_from_edg_list : function to find
% belonging triangle indices of the edge list.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2020-2025.
%
%
%%% Input arguments
%
%       [|  |  | ]
% - T = [i1 i2 i3], positive integer matrix double, the triangulation, size(T) = [nb_triangles,3]. Mandatory.
%       [|  |  | ]
%
%              [ | | ]
% - edg_list = [e1 e2], positive integer matrix double, the edge list, size(edg_list) = [nb_edg,2], with nb_edg the number of edges. Mandatory.
%              [ | | ]
%
%
%%% Output argument
%
% - tgl_id_list : positive integer row vector double, the triangle indices list, size(tgl_id_list) = [nb_triangles,3].


%% Body
% tic;
tgl_id_list = cellfun(@(r) find(sum(bitor(T == r(1,1),T == r(1,2)),2) == 2)',num2cell(edg_list,2),'un',0);
% fprintf('find_triangle_indices_from_edg_list request executed in %d seconds.\n',toc);


end % find_triangle_indices_from_edg_list