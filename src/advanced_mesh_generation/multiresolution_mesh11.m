function [V,T] = multiresolution_mesh11(V) % lvl % option surface  / volume
% multiresolution_mesh : function to mesh one given point set with triangular faces.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2020-2025.


% % Notes
%
% - split_edge ne peut pas être vectorisé comme splitrianglein3 car edge = élément commun à deux triangles;


% 0 Remove duplicated vertices
coeff = 1;
epsilon = coeff*eps;
V = uniquetol(V,epsilon,'ByRows',true);

% I Convex hull of the given point set inward oriented normals
T = convhull(V);
T = fliplr(T);


% I Nécessité d'un split_edge global, avec hystérésis
% Centré sur les edges, non sur les triangles
% -> Tout le problème vient de là (95%)
% -> une méthode utilisant build_triangulation_from_edg_list ?
%
% II + limites du fix_empty_octahedra -> peut créer des formes dégénérées...
%
% III Traiter les bord set les triangles
% par ordre de longueur / surface décroissante


new_vid = -1; % invalid but non empty vertex id to start with


% Main loop
% while ~isempty(new_vid)
for s = 1:1
                
    new_vid = [];
    
    N = face_normals(V,T);
    E = query_edg_list(T,'sorted');
    
    v_ban_id = unique(T(:))';    
    T = refine_edges(V, T, N, E, v_ban_id);
                    
    T = fliplr(T); % flip normals orientation
    N = face_normals(V,T);
    E = query_edg_list(T,'sorted');    
    [T,N] = fix_empty_octahedra(V,T,N,E,epsilon);    
    T = fliplr(T); % flip back inward triangles orientation
                    
    v_ban_id = unique(T(:))';
    T = refine_triangles(V,T,N,v_ban_id);
                    
    T = fliplr(T); % flip back inward normals orientation
    N = face_normals(V,T);
    E = query_edg_list(T,'sorted');    
    T = fix_empty_octahedra(V,T,N,E,epsilon);        
    T = fliplr(T); % flip back inward triangles orientation
   
end

T = fliplr(T); % flip back outward triangles orientation

% Remove duplicated triangles
T = remove_duplicated_triangles(T);


end % multiresolution_mesh


function T = refine_edges(V, T, N, E, v_ban_id)


N_e = edge_normals(T,N,E);
new_vk_id = [];


for e = 1:size(E,1)
    
    % Find closest point to the edge middle in
    % the opposite direction of edge normals
    
    I = mean(V(E(e,:),:),1); % edge middle
    dst = sqrt(sum((V-I).^2,2));
    cross_prod = cross(V-I,repmat(N_e(e,:)+I,[size(V,1),1]),2);
    % dot_prod = (V-I) * (N_e(e,:)+I)';
    criterion = sqrt(sum(cross_prod.^2,2)) .* dst;
    criterion(v_ban_id) = Inf;
    % criterion(dot_prod > -epsilon) = Inf;
    [~,new_vid] = min(criterion);
    
    if new_vid
        
        new_vid = new_vid(1,1);
        
        % Also if this point doesn't create
        % - Exterior / orphan point;
        % - Self intersecting triangles
        if isempty(find(v_ban_id == new_vid, 1))
            
            if isempty(find(new_vk_id == new_vid, 1))
                
                new_vk_id = cat(2,new_vk_id,new_vid);
                [~,T] = split_edge(V,T,E(e,:),'inset',new_vid); % ~ because no need of V which stays unchanged (no vertex addition)
                
            end
            
            v_ban_id = cat(2,v_ban_id,new_vid);
            
        end
        
    end
    
end


end


function T = refine_triangles(V, T, N, v_ban_id)


new_vk_id = [];

for tid = 1:size(T,1)
    
    % Find closest point to the face isobarycentre in
    % the opposite direction of face normals
    % G = mean(V(T(tid,:),:),1); % face/triangle isobarycentre
    
    [~,G] = triangle_circumcircle(V(T(tid,1),:)',V(T(tid,2),:)',V(T(tid,3),:)',3);
    dst = sqrt(sum((V-G').^2,2));
    cross_prod = cross(V-G',repmat(N(tid,:)+G',[size(V,1),1]),2);
    % dot_prod = (V-G') * (N(tid,:)+G)';
    criterion = sqrt(sum(cross_prod.^2,2)) .* dst;
    criterion(v_ban_id) = Inf;
    % criterion(dot_prod < epsilon) = Inf;
    [~,new_vid] = min(criterion);
    
    if new_vid
        
        new_vid = new_vid(1,1);
        
    end
    
    % Also if this point doesn't create
    % - Exterior / orphan point;
    % - Self intersecting triangles
    if isempty(find(v_ban_id == new_vid, 1))
        
        if isempty(find(new_vk_id == new_vid, 1))
            
            new_vk_id = cat(2,new_vk_id,new_vid);
            
        end
        
        v_ban_id = cat(2,v_ban_id,new_vid);
        
    end
    
end


if ~isempty(new_vk_id) % & numel(new_vk_id) ==size(T,1) % cas où pas de point trouvé pour certain triangles ?
    
    [~,T] = splitrianglein3(V,T,1:size(T,1),'insert',new_vk_id);
    
end


end % refine_triangles


function [T, N] = fix_empty_octahedra(V, T, N, E, epsilon)


eid = 1;

while eid < 1 + size(E,1)
    
    tgl_pair_id = cell2mat(find_triangle_indices_from_edg_list(T,E(eid,:)));
    
    if numel(tgl_pair_id) == 2
        
        isconvex = ~detect_concavity(V,T,N,tgl_pair_id,epsilon);
        
        if isconvex
            
            tgl_pair2 = find_tgl_tetra_id(tgl_pair_id,T);
            tetra = cat(1,T(tgl_pair_id,:),tgl_pair2);
            in_vtx_set_id = find(isin3Dconvexset(V,tetra,V,epsilon),1);
            
            if isempty(in_vtx_set_id)
                
                [T,N,E] = flip_two_ngb_triangles(tgl_pair_id,T,V,N,E);
                
            end
            
        end
        
    end
    
    eid = eid + 1;
    
end


end % fix_empty_octahedra


function isconcave = detect_concavity(V, T, N, tgl_pair_id, epsilon)
% detect_concavity : function to detect
% concave triangle pair configurations.


i1 = tgl_pair_id(1);
i2 = tgl_pair_id(2);

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


function [tgl_pair2] = find_tgl_tetra_id(tgl_pair_id, T)
% find_tgl_tetra_id


i1 = tgl_pair_id(1);
i2 = tgl_pair_id(2);

T1 = T(i1,:);
T2 = T(i2,:);

cmn_edg = intersect(T1,T2);
new_cmn_edg = setdiff(union(T1,T2),cmn_edg);

start_id1 = find(ismember(T1,new_cmn_edg(1)),1,'first');
T1c = circshift(T1,-1);

if ~isempty(start_id1)
    T3 = cat(2,new_cmn_edg(1),T1c(start_id1),new_cmn_edg(2));
else
    start_id1 = find(ismember(T1,new_cmn_edg(2)),1,'first');
    T3 = cat(2,new_cmn_edg(2),T1c(start_id1),new_cmn_edg(1));
end


start_id2 = find(ismember(T2,new_cmn_edg(2)),1,'first');
T2c = circshift(T2,-1);

if ~isempty(start_id2)
    T4 = cat(2,new_cmn_edg(2),T2c(start_id2),new_cmn_edg(1));
else
    start_id2 = find(ismember(T2,new_cmn_edg(1)),1,'first');
    T4 = cat(2,new_cmn_edg(1),T2c(start_id2),new_cmn_edg(2));
end

tgl_pair2 = cat(1,fliplr(T3),fliplr(T4));


end % find_tgl_tetra_id


% TODO : use flip_edge -> voir enveloppe convexe
function [T, N, edg_list] = flip_two_ngb_triangles(tgl_pair_id, T, V, N, edg_list)
% flip_two_ngb_triangles : function to flip
% two triangles sharing one common edge.


T1 = T(tgl_pair_id(1),:);
T2 = T(tgl_pair_id(2),:);

cmn_edg = intersect(T1,T2);                  % = new opposit vertices
new_cmn_edg = setdiff(union(T1,T2),cmn_edg); % = current opposit vertices


start_id1 = find(ismember(T1,new_cmn_edg(1)),1,'first');
T1c = circshift(T1,-1);

if ~isempty(start_id1)
    Ta = cat(2,new_cmn_edg(1),T1c(start_id1),new_cmn_edg(2));
else
    start_id1 = find(ismember(T1,new_cmn_edg(2)),1,'first');
    Ta = cat(2,new_cmn_edg(2),T1c(start_id1),new_cmn_edg(1));
end


start_id2 = find(ismember(T2,new_cmn_edg(2)),1,'first');
T2c = circshift(T2,-1);

if ~isempty(start_id2)
    Tb = cat(2,new_cmn_edg(2),T2c(start_id2),new_cmn_edg(1));
else
    start_id2 = find(ismember(T2,new_cmn_edg(1)),1,'first');
    Tb = cat(2,new_cmn_edg(1),T2c(start_id2),new_cmn_edg(2));
end


% TODO : substituer (conserver le même index de ligne)

% Add 2 triangles and their face normals
T = add_triangles(Ta,T);
T = add_triangles(Tb,T);

% Remove 2 triangles and their face normals
T(tgl_pair_id,:) = [];
new_face_normals = face_normals(V,T(end-1:end,:),'norm');
N = cat(1,N,new_face_normals);
N(tgl_pair_id,:) = [];
edg_list = cat(1,edg_list,sort(new_cmn_edg));
edg_list(all(bsxfun(@eq,edg_list,sort(cmn_edg)),2),:) = [];


end % flip_two_ngb_triangles


% % dst = sqrt(sum((V-I).^2,2)); % rajouter "dans la direction opposée à sa normale" -> produit scalaire...


% dot_prod = (V - I) * N_e(e,:)';
% dot_prod(v_ban_id) = -Inf;


% dot_prod = (V-I) * N_e(e,:)';
% cross_prod = cross(V-I,repmat(N_e(e,:),[size(V,1),1]),2);
% criterion = dot_prod ./ cross_prod;
% criterion(v_ban_id) = -Inf;
% [~,new_vid] = max(criterion);

% Ou sinon min produit scalaire et min produit vectoriel (somme ou produit ?)