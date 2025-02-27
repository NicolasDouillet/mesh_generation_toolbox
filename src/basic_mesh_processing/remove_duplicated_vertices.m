function [V_out, T_out] = remove_duplicated_vertices(V_in, T_in, tol)
%% remove_duplicated_vertices : function to remove duplicated
% vertices from the point set.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2020-2025.
%
%
%%% Input arguments
%
%          [ |    |    |  ]
% - V_in = [X_in Y_in Z_in], real matrix double, the input point set, size(V_in) = [nb_input_vertices,3]. Mandatory.
%          [ |    |    |  ]
%
%          [  |     |     |  ]
% - T_in = [i1_in i2_in i3_in], positive integer matrix double, the input triangulation, size(T_in) = [nb_input_triangles,3]. Mandatory.
%          [  |     |     |  ]
%
% - tol : real scalar double. tol = f(eps), the tolerance to error. Mandatory.
%
%
%%% Output arguments
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
% tic;
[V_out,~,n] = uniquetol(V_in,tol,'ByRows',true);
T_out = n(T_in);
% fprintf('%d duplicated vertices removed in %d seconds.\n',size(V_in,1)-size(V_out,1),toc);


end % remove_duplicated_vertices