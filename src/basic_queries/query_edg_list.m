function edg_list = query_edg_list(T, mode)
%% query_edg_list : function to query the edges list
% corresponding to the triangulation T.
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
% - mode : character string in the set {'sorted','raw','SORTED','RAW'}, the
%          variable deciding wether to sort or not the edges in ascending
%          order. Case insensitive. Optional.
%
%
%%% Output argument
%
%              [ | | ]
%              [i1 i2]
% - edg_list = [i2 i3], positive integer matrix double, the edges list, size(edg_list) =  [nb_edg,2]
%              [i3 i1]
%              [ | | ]
%
%              with nb_edg the number of edges.


%% Body
% tic;
L = cat(2,T,T(:,1)); % loop
R = repelem(L,1,[1 2 2 1]); % replicated
face_nb_vtx = size(T,2);

edg_list = unique(cell2mat(cellfun(@(x) reshape(x,[2,face_nb_vtx])',num2cell(R,2),'un',0)),'rows');

if nargin  > 1 && strcmpi(mode,'sorted')
    
    edg_list = unique(sort(edg_list,2),'rows');
    
    % elseif nargin  < 2 || strcmpi(mode,'raw')
    
    % do nothing
    
    % else
    
    % do nothing
    
end

% fprintf('query_edg_list request executed in %d seconds.\n',toc);


end % query_edg_list