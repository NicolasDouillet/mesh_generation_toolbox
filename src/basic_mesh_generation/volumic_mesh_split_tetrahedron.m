function [V, T] = volumic_mesh_split_tetrahedron(nb_it, V1, V2, V3, V4, sampling_mode)
%% volumic_mesh_split_tetrahedron : function to build and split a volumic mesh of a given
% tetrahedron. Performs either middle-edge or face-barycentric subdivisions.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2025.
%
%
%%% Input arguments
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
% - sampling_mode : character string in the set {'barycentric','BARYCENTRIC','edge_subdivision'*,'EDGE_SUBDIVISION'}, the sampling_mode. Case insensitive. Optional.
%
%
%%% Output arguments
%
%       [|  |  | ]
% - V = [Vx Vy Vz], real matrix double, the output point set, size(V) = [nb_output_vertices,3],
%       [|  |  | ]
%
%       [|  |  | ]
% - T = [t1 t2 t3], positive integer matrix double, the triangulation, size(T) = [nb_triangles,3].
%       [|  |  | ]


%% Input parsing
if nargin < 6
    
    sampling_mode = 'barycentric';
    
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
        
        if strcmpi(sampling_mode,'barycentric')
            
            [V,T,new_tetra] = split_tetra_barycentre(V,T,tetra_list(k,:));
        
        elseif strcmpi(sampling_mode,'edge_subdivision')
            
            [V,T,new_tetra] = split_tetra_edges(V,T,tetra_list(k,:));
            
        end
        
        T = remove_duplicated_triangles(T);
        tetra_list = cat(1,tetra_list,new_tetra);
        
    end
    
    tetra_start_idx = nb_tetra + 1;
    
end


end % volumic_mesh_split_tetrahedron


%% split_tetra_barycentre subfunction
function [V, T, new_tetra] = split_tetra_barycentre(V, T, tetra_list_k) 
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2024-2025.


% One new vertex
new_vtx = mean(V(tetra_list_k,:),1); 
V = cat(1,V,new_vtx);
new_vtx_idx = size(V,1);

% 6 new triangles
new_tgl = combnk([tetra_list_k,new_vtx_idx],3);
T = cat(1,T,new_tgl);

% 4 new tetrahedra
new_tetra = cat(2,combnk(tetra_list_k,3),new_vtx_idx*ones(4,1)); 


end % split_tetra_barycentre


%% split_tetra_edges subfunction
function [V, T, new_tetra] = split_tetra_edges(V, T, tetra_list_k) 
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2024-2025.


% 7 new vertices
% Middle edge ones
new_vtx1 =  mean(V([tetra_list_k(1),tetra_list_k(2)],:),1);
new_vtx2 =  mean(V([tetra_list_k(1),tetra_list_k(3)],:),1);
new_vtx3 =  mean(V([tetra_list_k(1),tetra_list_k(4)],:),1);
new_vtx4 =  mean(V([tetra_list_k(2),tetra_list_k(3)],:),1);
new_vtx5 =  mean(V([tetra_list_k(2),tetra_list_k(4)],:),1);
new_vtx6 =  mean(V([tetra_list_k(3),tetra_list_k(4)],:),1);

new_vtx7 =  mean(V(tetra_list_k,:),1); % central one

V = cat(1,V,new_vtx1,new_vtx2,new_vtx3,new_vtx4,new_vtx5,new_vtx6,new_vtx7);
new_vtx_idx = size(V,1)-6:size(V,1);

% 12 new tetrahedra
new_tetra = [tetra_list_k(1) new_vtx_idx(1) new_vtx_idx(2) new_vtx_idx(3);...
             tetra_list_k(2) new_vtx_idx(1) new_vtx_idx(4) new_vtx_idx(5);...
             tetra_list_k(3) new_vtx_idx(2) new_vtx_idx(4) new_vtx_idx(6);...
             tetra_list_k(4) new_vtx_idx(3) new_vtx_idx(5) new_vtx_idx(6);...
             
             new_vtx_idx(1) new_vtx_idx(2) new_vtx_idx(4) new_vtx_idx(7);...
             new_vtx_idx(4) new_vtx_idx(5) new_vtx_idx(6) new_vtx_idx(7);...                                                       
             new_vtx_idx(1) new_vtx_idx(3) new_vtx_idx(5) new_vtx_idx(7);...
             new_vtx_idx(2) new_vtx_idx(3) new_vtx_idx(6) new_vtx_idx(7);...   
             
             new_vtx_idx(1) new_vtx_idx(2) new_vtx_idx(3) new_vtx_idx(7);...
             new_vtx_idx(1) new_vtx_idx(4) new_vtx_idx(5) new_vtx_idx(7);...                                                       
             new_vtx_idx(2) new_vtx_idx(4) new_vtx_idx(6) new_vtx_idx(7);...
             new_vtx_idx(3) new_vtx_idx(5) new_vtx_idx(6) new_vtx_idx(7);...
             ];

% 48 new triangles
new_tgl = zeros(48,3);

for n = 1:size(new_tetra,1)

    new_tgl(4*n-3:4*n,:) = combnk(new_tetra(n,:),3);
    
end

T = cat(1,T,new_tgl);


end % split_tetra_edges