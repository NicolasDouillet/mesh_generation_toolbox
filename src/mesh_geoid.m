function [V, T] = mesh_geoid(id, nb_it, sampling_mode)
% mesh_geoid : function to mesh a geoid based one platonic solid
% (all except dodecahedron) iterative projections on the unit sphere,
% with two different sampling modes available.
%
% Author & support : nicolas.douillet (at) free.fr, 2021-2023.
%
%
%%% Input arguments
%
% - id :              positive integer scalar double, the basis polyhedron (platonic solid) id.
%
% - sampling_mode :   character string in the set : {'edge','face'}. Case insensitive.
%
% - nb_it :           positive integer scalar double, either the number of sub edges to subdivide the original edge in,
%                     when sampling_mode = 'edge') or the number of iterations to perform,
%                     when sampling_mode = 'face'.
%
%
%%% Output arguments
%
%       [| | |]
% - V = [X Y Z], real matrix double, the resulting point set, size(V) = [nb_vertices,3].
%       [| | |]
%
%       [ |  |  |]
% - T = [i1 i2 i3], positive integer matrix double, the resulting triangulation, size(T) = [nb_triangles,3].
%       [ |  |  |]
%
%
%%% About / other informations
%
% Geoid is centered on the origin, [0 0 0].
% Triangles / normals are coherently oriented and facing outward.


epsilon = 1e3*eps;
[V,T] = platonic_solids(id,1,'triangle');


if nargin < 3
   
    sampling_mode = 'edge';
    
end


switch id
    
    case 1 % tetrahedron
        
        begin_nb_faces = 4;
        
    case 2 % triangulated cube
        
        begin_nb_faces = 12;
        
    case 3 % octahedron
        
        begin_nb_faces = 8;
        
    case 4 % icosahedron
        
        begin_nb_faces = 20;
        
    case 5 % dodecahedron
        
        begin_nb_faces = 36;
        
%     otherwise
%         
%         error('id must be in the range |[1 ; 5]|.');
        
end

        
if strcmpi(sampling_mode,'edge')
    
    % Oversample triangles by creating new vertices ; link vertices to create new triangles
    nt = nb_it^2;               % nb new triangles
    nv = (nb_it+1)*(nb_it+2)/2; % nb new vertices (1:nb_it) sum
    
    T_new = zeros(begin_nb_faces*nt,3);
    V_new = zeros(0.5*begin_nb_faces*(nb_it+1)*(nb_it+2),3);
    
    for j = 1:size(T,1)
        
        [new_sub_V, new_sub_T] = mesh_triangle(V(T(j,1),:)',V(T(j,2),:)',V(T(j,3),:)',nb_it);
        new_sub_T = new_sub_T + (j-1)*nv; % update triangle indices
        
        T_new((j-1)*nt+1:j*nt,:) = new_sub_T;
        V_new((j-1)*nv+1:j*nv,:) = new_sub_V;
        
    end
    
    T = T_new;
    V = V_new;
    
    % Vertices normalization / projection on the sphere surface
    V = V ./ sqrt(sum(V.^2,2));        
    
elseif strcmpi(sampling_mode,'face')
    
    iteration = 0; % to start with
    N = compute_face_normals(V,T,'norm');
    
    while iteration < nb_it
        
        tgl_idx = 1:begin_nb_faces*3^iteration; % concavity
        
        [V,T,N] = grow_nxt_lvl_tetrahedra(V,T,N,tgl_idx);
        edg_list = query_edges_list(T,'sorted');
        i = 1;
        
        while i < 1 + size(edg_list,1)
            
            tgl_pair_idx = cell2mat(find_triangle_indices_from_edges_list(T,edg_list(i,:)));
            isconcave = detect_concavity(V,T,N,tgl_pair_idx,epsilon);
            
            if isconcave
                
                [T,N,edg_list] = flip_two_ngb_triangles(tgl_pair_idx,T,V,N,edg_list);
                
            else
                
                i = i + 1;
                
            end
            
        end
        
        iteration = iteration + 1;
        
    end
    
else
    
    error('Unknown sampling mode.');
    
end


% Remove duplicated vertices
[V,~,n] = uniquetol(V,epsilon,'ByRows',true);
T = n(T);

% Remove duplicated triangles
T_sort = sort(T,2);
[~,idx,~] = unique(T_sort,'rows','stable');
T = T(idx,:);

% Reorient normal outward
T = fliplr(T);


end % mesh_geoid


% grow_nxt_lvl_tetrahedra subfunction
function [V, T, N] = grow_nxt_lvl_tetrahedra(V, T, N, tgl_idx)
% grow_nxt_lvl_tetrahedra : function to create the three new
% vertices and link them to the 4*3^nb_it new triangles
% to create the next triangulation level and erase
% the 4*three^(nb_it-1) previous triangles.
%
% Author & support : nicolas.douillet (at) free.fr, 2021-2023.


for idx = tgl_idx
    
    new_vtx = mean(V(T(idx,:),:),1);
    new_vtx = new_vtx ./ sqrt(sum(new_vtx.^2,2));
    V = cat(1,V,new_vtx);
    
    new_tgl1 = cat(2,T(idx,1:2),size(V,1));
    new_tgl2 = cat(2,T(idx,2:3),size(V,1));
    new_tgl3 = cat(2,T(idx,3),T(idx,1),size(V,1));
    
    % Add 3 new triangles and face normals
    T = cat(1,T,new_tgl1,new_tgl2,new_tgl3);
    new_face_normals = compute_face_normals(V,T(end-2:end,:),'norm');
    N = cat(1,N,new_face_normals);
    
end

% Remove one triangle and its normal
T(tgl_idx,:) = [];
N(tgl_idx,:) = [];


end % grow_nxt_lvl_tetrahedra


% detect_concavity subfunction
function [isconcave] = detect_concavity(V, T, N, tgl_pair_idx, epsilon)
% detect_concavity : function to detect concave triangle pair configurations.
%
% Author & support : nicolas.douillet (at) free.fr, 2021-2023.


i1 = tgl_pair_idx(1);
i2 = tgl_pair_idx(2);

T1 = T(i1,:);
T2 = T(i2,:);

cmn_edg = intersect(T1,T2);
cross_edg = setdiff(union(T1,T2),cmn_edg);

H1 = mean(V(cmn_edg,:),1);
H2 = mean(V(cross_edg,:),1);

n1 = N(i1,:);
n2 = N(i2,:);

isconcave = sign(dot(n1+n2,H2-H1,2).*(abs(dot(n1+n2,H2-H1,2)) > epsilon ) ) > 0;


end % detect_concavity


% flip_two_ngb_triangles subfunction
function [T, N, edg_list] = flip_two_ngb_triangles(tgl_pair_idx, T, V, N, edg_list)
% flip_two_ngb_triangles : function to flip two triangles sharing one common edge.
%
% Author & support : nicolas.douillet (at) free.fr, 2021-2023.


nb_vtx = size(V,1);
T1 = T(tgl_pair_idx(1),:);
T2 = T(tgl_pair_idx(2),:);

cmn_edg = intersect(T1,T2);                  % = new opposit vertices
new_cmn_edg = setdiff(union(T1,T2),cmn_edg); % = current opposit vertices


start_idx1 = find(ismember(T1,new_cmn_edg(1)),1,'first');
T1c = circshift(T1,-1);

if ~isempty(start_idx1)
    Ta = cat(2,new_cmn_edg(1),T1c(start_idx1),new_cmn_edg(2));
else
    start_idx1 = find(ismember(T1,new_cmn_edg(2)),1,'first');
    Ta = cat(2,new_cmn_edg(2),T1c(start_idx1),new_cmn_edg(1));
end


start_idx2 = find(ismember(T2,new_cmn_edg(2)),1,'first');
T2c = circshift(T2,-1);

if ~isempty(start_idx2)
    Tb = cat(2,new_cmn_edg(2),T2c(start_idx2),new_cmn_edg(1));
else
    start_idx2 = find(ismember(T2,new_cmn_edg(1)),1,'first');
    Tb = cat(2,new_cmn_edg(1),T2c(start_idx2),new_cmn_edg(2));
end

% Add 2 triangles and their face normals
T = add_triangles(Ta,T,nb_vtx);
T = add_triangles(Tb,T,nb_vtx);

% Remove 2 triangles and their face normals
T(tgl_pair_idx,:) = [];
new_face_normals = compute_face_normals(V,T(end-1:end,:),'norm');
N = cat(1,N,new_face_normals);
N(tgl_pair_idx,:) = [];
edg_list = cat(1,edg_list,sort(new_cmn_edg));
edg_list(all(bsxfun(@eq,edg_list,sort(cmn_edg)),2),:) = [];


end % flip_two_ngb_triangles