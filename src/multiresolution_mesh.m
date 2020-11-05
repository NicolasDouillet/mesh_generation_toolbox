function [T] = multiresolution_mesh(V, lvl)


% surf_type = 'closed'; % 'opened'

% principle : at each iteration look for and choose min(angle(n,GhVi) x ||GhVi||)
% ou juste min distance sinon
% -> dig_tetrahedron


epsilon = 1e8*eps; % floating point tolerance error

T = convhull(V);
nb_t = size(T,1);
Vc = V(unique(T(:)),:);
N = compute_face_normals(V,T,'norm');
G = mean(Vc,1);
Gt = cell2mat(cellfun(@(r) mean(V(r,:),1),num2cell(T,2),'un',0));

% Orient normals inward
orientation = sign(dot(N,Gt-repmat(G,[nb_t,1]),2));

if ~isequal(orientation,-ones(nb_t,1)) && ~isequal(orientation,zeros(nb_t,1))
    N = N.*orientation;
    T(orientation > 0,:) = fliplr(T(orientation > 0,:));
end


nb_it = 0;


while nb_it < lvl
    
    curr_tgl_idx = 1;
    
    while curr_tgl_idx < 3 % 1 + size(T,1)
        
        [T,N,new_vtx_idx] = dig_tetrahedron(V,T,N,curr_tgl_idx,epsilon);
        
        %         if new_vtx_idx % effective dig with new triangles
        %
        %             edg_list = query_edges_list(T,'sorted');
        %             i = 1;
        %
        %             while i < 1 + size(edg_list,1)
        %
        %                 tgl_pair_idx = cell2mat(find_triangle_indices_from_edges_list(T,edg_list(i,:)));
        %
        %                 if numel(tgl_pair_idx) == 2
        %
        %                     isconcave = detect_concavity(V,T,N,tgl_pair_idx,epsilon); % TODO : flat cases management
        %                     
        %                     
        %                     % remove_triangles_with_lonely_edge
        %                     
        %
        %
        %                     if isconcave % = convex with inverse normals
        %
        %                         tgl_pair2 = find_tgl_tetra_idx(tgl_pair_idx,T);
        %                         tetra = cat(1,T(tgl_pair_idx,:),tgl_pair2);
        %                         in_vtx_set_idx = find(isin3Dconvexset(V,tetra,V,epsilon)); % TODO : flat cases management
        %
        %                         if isempty(in_vtx_set_idx)
        %
        %                             disp('flip_two_ngb_triangles');
        %                             [T,N,edg_list] = flip_two_ngb_triangles(tgl_pair_idx,T,V,N,edg_list);
        %
        %                         end
        %
        %                         % else
        %
        %
        %                     end
        %
        %                 else
        %
        %                     break;
        %
        %                 end
        %
        %                 i = i + 1;
        %
        %             end
        %
        %             % curr_tgl_idx = curr_tgl_idx - 1;
        %
        %         end
        
        curr_tgl_idx = curr_tgl_idx + 1;
        
    end
    
    nb_it = nb_it + 1;
    
end


% % In case of opened surface
% remove duplicated triangles
% remove self intersecting triangles ?


end % multiresolution_mesh


function [T, N, new_vtx_idx] = dig_tetrahedron(V, T, N, tgl_idx, epsilon)
% dig_tetrahedron : function to find one
% new vertex belonging to the triangulation
% for one given triangle, to create the three
% newborn triangles and erase the previous one.


vect_angle = @(cross_prod,dot_prod) atan2(sqrt(sum(cross_prod.^2,2)),dot_prod);

nb_vtx = size(V,1);


% New criterion = min(d(Gh,H)) * angle(GhVi,N)) * dot_prod(GhVi,Ni))


Gh = mean(V(T(tgl_idx,:),:),1);
GhVi = V - repmat(Gh,[nb_vtx,1]);
Ni = repmat(N(tgl_idx,:),[nb_vtx,1]);
v_angle = vect_angle(cross(Ni,GhVi,2),dot(Ni,GhVi,2));
[d2Hi,Hi] = point_to_plane_distance(V,N(tgl_idx,:),Gh);
dst = sqrt(sum((Hi-repmat(Gh,[nb_vtx,1])).^2,2));


criterion = dst .* v_angle .* d2Hi;
criterion(criterion < epsilon) = Inf; % /_!_\ ne pas supprimer de valeur pour garder les bons indices /_!_\

[~,mink_idx] = mink(criterion,size(criterion,1));
f = setdiff(mink_idx,unique(T(:))); % avoid non manifold vertex creation
f = f(1,1);

new_tgl1 = cat(2,T(tgl_idx,1:2),f);
new_tgl2 = cat(2,T(tgl_idx,2:3),f);
new_tgl3 = cat(2,T(tgl_idx,3),T(tgl_idx,1),f);


% Forbid duplicated triangles
% Forbid flat triangles


% Add 3 new triangles and face normals
T = cat(1,T,new_tgl1,new_tgl2,new_tgl3);
new_face_normals = compute_face_normals(V,T(end-2:end,:),'norm');
N = cat(1,N,new_face_normals);

% Remove one triangle and its normal
T(tgl_idx,:) = [];
N(tgl_idx,:) = [];
new_vtx_idx = f;


end % dig_tetrahedron


function [tgl_pair2] = find_tgl_tetra_idx(tgl_pair_idx, T)
% find_tgl_tetra_idx


i1 = tgl_pair_idx(1);
i2 = tgl_pair_idx(2);

T1 = T(i1,:);
T2 = T(i2,:);

cmn_edg = intersect(T1,T2);
new_cmn_edg = setdiff(union(T1,T2),cmn_edg);

start_idx1 = find(ismember(T1,new_cmn_edg(1)),1,'first');
T1c = circshift(T1,-1);

if ~isempty(start_idx1)
    T3 = cat(2,new_cmn_edg(1),T1c(start_idx1),new_cmn_edg(2));
else
    start_idx1 = find(ismember(T1,new_cmn_edg(2)),1,'first');
    T3 = cat(2,new_cmn_edg(2),T1c(start_idx1),new_cmn_edg(1));
end


start_idx2 = find(ismember(T2,new_cmn_edg(2)),1,'first');
T2c = circshift(T2,-1);

if ~isempty(start_idx2)
    T4 = cat(2,new_cmn_edg(2),T2c(start_idx2),new_cmn_edg(1));
else
    start_idx2 = find(ismember(T2,new_cmn_edg(1)),1,'first');
    T4 = cat(2,new_cmn_edg(1),T2c(start_idx2),new_cmn_edg(2));
end

tgl_pair2 = cat(1,fliplr(T3),fliplr(T4));


end % find_tgl_tetra_idx


function [isconcave] = detect_concavity(V, T, N, tgl_pair_idx, epsilon)
% detect_concavity : function to detect
% concave triangle pair configurations.


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

isconcave = sign(dot(n1+n2,H2-H1,2)) > epsilon;


end % detect_concavity


function [T, N, edg_list] = flip_two_ngb_triangles(tgl_pair_idx, T, V, N, edg_list)
% flip_two_ngb_triangles : function to flip
% two triangles sharing one common edge.


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