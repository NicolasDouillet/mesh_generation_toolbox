function [V_out, T_out] = remove_vertices(V_set, V_in, T_in, mode)
%% remove_vertices : function to remove vertices from the vertex set.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2020-2025.
%
%
%%% Input arguments
%
%           [| | |]
% - V_set = [X Y Z], real matrix double, the vertex set to remove, size(V_set) = [nb_new_vertices,3],
%           [| | |]
%                    or positive integer row vector double, the index list of the vertices to remove. Mandatory.
%
%          [ |    |    |  ]
% - V_in = [X_in Y_in Z_in], real matrix double, the input point set, size(V_in) = [nb_input_vertices,3]. Mandatory.
%          [ |    |    |  ]
%
%          [  |     |     |  ]
% - T_in = [i1_in i2_in i3_in], positive integer matrix double, the input triangulation, size(T_in) = [nb_input_triangles,3]. Mandatory.
%          |  |     |     |  ]
%
% - mode : character string in the set {'index','explicit','INDEX','EXPLICIT'}. Explicit mode corresponds to the
%          case where V_set is made of vertices [X Y Z] coordinates. Case insensitive. Optional.
%
%
%%% Output arguments
%
%           [  |     |     |  ]
% - V_out = [X_out Y_out Z_out], real matrix double, the output point set, size(V_out) = [nb_output_vertices,3],
%           [  |     |     |  ]
%
%           with nb_output_vertices = nb_input_vertices + nb_new_vertices.
%
%           [  |      |      |   ]
% - T_out = [i1_out i2_out i3_out], positive integer matrix double, the output triangulation, size(T_out) = [nb_output_triangles,3].
%           [  |      |      |   ]


% tic;

%% Input parsing
if nargin  < 4
    
    mode = 'index';
    
else
    
    if ~(strcmpi(mode,'index') && ismember(1,size(V_set))) && ~(strcmpi(mode,'explicit') && size(V_set,2) == 3)
           
        error('mode value must be either set to ''index'' with V_set a one row/column indices vector or to ''explicit'' with size(V_set,2) = 3.');
       
    end
    
end


%% Body
V_out = V_in;
T_out = T_in;

if strcmpi(mode,'index') && ismember(1,size(V_set)) || nargin < 4
       
    vtx_id_2_rm = V_set;                   
    
elseif strcmpi(mode,'explicit') && size(V_set,2) == 3
       
    vtx_id_2_rm = ismember(V_in,V_set,'rows');
    
end

% I Suppress vertices & triangles
V_out(vtx_id_2_rm,:) = [];
tgl_id_2_rm = find_triangles_from_vertex_list(T_in,vtx_id_2_rm);
T_out(tgl_id_2_rm,:) = [];

% II Reindex triangulation
keep_id = setdiff(1:size(V_in,1),vtx_id_2_rm); % what V_out indices are relatively to V_in ones
M = containers.Map(keep_id,1:size(V_out,1));
T_suppr_id = isKey(M,num2cell(T_out));
T_out(T_suppr_id) = cell2mat(values(M,num2cell(T_out(T_suppr_id))));

% fprintf('%d vertices removed in %d seconds.\n',nnz(vtx_id_2_rm),toc);


end % remove_vertices