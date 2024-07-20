function [V_out, new_id] = add_vertices(V_set, V_in)
%% add_vertices : function to add some new vertices to the vertex set.
%
% Author : nicolas.douillet9 (at) gmail.com, 2020-2024.
%
%
% Input arguments
%
%           [| | |]
% - V_set = [X Y Z], real matrix double, the vertex set to add, size(V_set) = [nb_new_vertices,3].
%           [| | |]
%
%          [ |    |    |  ]
% - V_in = [X_in Y_in Z_in], real matrix double, the input point set, size(V_in) = [nb_input_vertices,3].
%          [ |    |    |  ]
%
%
%
% Output arguments
%
%           [  |     |     |  ]
% - V_out = [X_out Y_out Z_out], real matrix double, the output point set, size(V_out) = [nb_output_vertices,3],
%           [  |     |     |  ]
%
%           with nb_output_vertices = nb_input_vertices + nb_new_vertices.
%
% - new_id : positive integer row vector double, added vertices indices, size(new_id) = [1,nb_new_vtx].


%% Body
% tic;
coeff = 1;
epsilon = coeff*eps;

if isreal(V_set)
    
    dpl_vtx_id = ismembertol(sort(V_set,2),sort(V_in,2),epsilon,'ByRows',true);
    
    if nnz(dpl_vtx_id)
        
        % Suppress duplicated vertices
        V_set = V_set(~dpl_vtx_id,:);
        % warning('One or more vertex from this set already exist in the current vertex set. Duplicated vertices have been ignored.\n');
        
    end
    
    V_out = cat(1,V_in,V_set);
    new_id = (1+size(V_in,1)):size(V_in,1)+size(V_set,1);
    
else
    
    error('Unable to perform vertex addition because V_set contains some invalid vertex coordinates. Vertex coordinates must be real numbers.');
    
end

% fprintf('%d vertices added in % d seconds.\n',size(V_set,1),toc);


end % add_vertices