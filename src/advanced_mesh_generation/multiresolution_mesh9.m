function [V, T_new] = multiresolution_mesh9(V, lvl)
% multiresolution_mesh : function to mesh one given point set with triangular faces.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2024-2025.


% -> "splitedge" (avec point différent du milieu) dans la direction de la normale orientée vers l'intérieure.
% Boucle sur les edges par ordre décroissant de longueur.
% Sans mise à jour à chaque itération, mais mise à jour edge list après chaque boucle.

% Nécessite de calculer toutes les edge normals orientée vers l'intérieur
% (calculées à partir des normales des faces normées ou non ?)

% Faces en dernier

% -> Méthode pour ne recalculer que partiellement les vertex et face normals
%
% - update_vertex_normals
% - update_face_normals



% Comment s'assurer de ne pas "traverser un point" ?



addpath('C:\Users\Nicolas\Desktop\TMW_contributions\mesh_processing_toolbox\src\basic_mesh_processing\'); % for splitrianglein3
addpath('C:\Users\Nicolas\Desktop\TMW_contributions\mesh_processing_toolbox\src\geometry_ressources\');   % for isin3Dtriangle


% I Convex hull with inward oriented normals
T = convhull(V);
T = fliplr(T);
N = face_normals(V,T,'norm');
T_cur = T;
T_new = T_cur; % special case of lvl = 0
N_cur = N;

coeff = 1;
epsilon = coeff*eps;


% II Loop on lvl number of steps
for k = 1:lvl % or max number of faces as a criterion
        
    G_i = cell2mat(cellfun(@(r) mean(V(r,:),1),num2cell(T_cur,2),'un',0));
    % V_ban_id = unique(T_cur(:)); % convex hull vertices are forbidden
           
    T_new = T_cur;
    N_new = N_cur;
    tgl_id2rm = [];
    n = 1;
    
     % III Loop on every triangle of the current mesh
    while n < 1 + size(T_cur,1)
        
        
        % IV Seek closest point to the triangle surface and which projects inside it
        % => should avoid self intersecting triangles
        %
        % -> point_to_plane_distance
        % -> isin3Dtriangle
        %
        % TODO : distance signée  = produit scalaire entre Ni et HiV) min positif ou nul                
        
        
        V_ban_id = T_new(n,:); % convex hull vertices are forbidden
        
        %         [d2Hi,Hi] = point_to_plane_distance(V,N_new(n,:),V(T_new(n,1),:));
        %         isin = zeros(size(V,1),1);
        %
        %         for q = 1:size(V,1)
        %
        %             isin(q,1) = isPointin3Dtriangle(V(T_new(n,1),:),V(T_new(n,2),:),V(T_new(n,3),:),Hi(q,:));
        %
        %         end
        %
        %         cd_vtx_id = setdiff(find(isin),T_new(n,:));
        %         % cd_vtx_sgn_dst2tgl = dot(repmat(N_new(n,:),[numel(cd_vtx_id),1]),V(cd_vtx_id,:)-Hi(cd_vtx_id,:));
        %
        %         % Distance non signée ici !
        %         [~,id] = min(d2Hi(cd_vtx_id,1));
        %
        %         if numel(id) > 1
        %
        %             id = id(1,1);
        %
        %         end
        %
        %         nrst_vtx_id = cd_vtx_id(id);
        
        
            dst2Gi = sqrt(sum((G_i(n,:)-V).^2,2));
            [~,vtx_id] = min(dst2Gi); % 2 in order to always be able to find a closest non banned vertex
            nrst_vtx_id = setdiff(vtx_id,V_ban_id);

            if numel(nrst_vtx_id) > 1
    
                nrst_vtx_id = nrst_vtx_id(1,1);
    
            end
        
        
        if ~isempty(nrst_vtx_id)
            
            V_ban_id = cat(2,V_ban_id,nrst_vtx_id); % to avoid non manifold vertices creation
            
            % New triangles and normals
            new_tgl1 = cat(2,T_new(n,1:2),nrst_vtx_id);
            new_tgl2 = cat(2,T_new(n,2:3),nrst_vtx_id);
            new_tgl3 = cat(2,T_new(n,3),T_new(n,1),nrst_vtx_id);
            
            new_N1 = cross(V(new_tgl1(2),:)-V(new_tgl1(1),:),V(new_tgl1(end),:)-V(new_tgl1(1),:),2);
            new_N1 = new_N1./sqrt(sum(new_N1.^2,2));
            
            new_N2 = cross(V(new_tgl2(2),:)-V(new_tgl2(1),:),V(new_tgl2(end),:)-V(new_tgl2(1),:),2);
            new_N2 = new_N2./sqrt(sum(new_N2.^2,2));
            
            new_N3 = cross(V(new_tgl3(2),:)-V(new_tgl3(1),:),V(new_tgl3(end),:)-V(new_tgl3(1),:),2);
            new_N3 = new_N3./sqrt(sum(new_N3.^2,2));
            
            % Update triangulation and normals
            T_set = cat(1,new_tgl1,new_tgl2,new_tgl3);
            T_new = cat(1,T_new,T_set);
            
            N_set = cat(1,new_N1,new_N2,new_N3);
            N_new = cat(1,N_new,N_set);
            
            tgl_id2rm = cat(1,tgl_id2rm,n);
            
        else
            
            % Nothing to do
            
        end
                
        n = n + 1;                        
        
    end
    
    
    % Remove previous triangles and normals
    T_new(tgl_id2rm,:) = [];
    T_cur = T_new;
    N_new(tgl_id2rm,:) = [];
    N_cur = N_new;
    
    % Detect concavity + flip empty triangle pairs
    % Or to perform after the creation of each new triangle ?
    E = query_edg_list(T_cur,'sorted');
    i = 1;
    
    while i < size(E,1)
        
        tgl_pair_id = cell2mat(find_triangle_indices_from_edg_list(T_cur,E(i,:)));
        
        if numel(tgl_pair_id) == 2
            
            isconcave = detect_concavity(V,T_cur,N_cur,tgl_pair_id,epsilon); 
            
            
            % -> Plus besoin lorsque l'algo sera "propre et optimal"
            if isconcave % = convex with inverse normals
                
                tgl_pair2 = find_tgl_tetra_id(tgl_pair_id,T_cur);
                tetra = cat(1,T_cur(tgl_pair_id,:),tgl_pair2);
                in_vtx_set_id = find(isin3Dconvexset(V,tetra,V,epsilon),1);
                
                if isempty(in_vtx_set_id)
                                        
                    [T_cur,N_cur,E] = flip_two_ngb_triangles(tgl_pair_id,T_cur,V,N_cur,E);
                    
                end                                
                
            end
            
        else
            
            break;
            
        end
        
        i = i + 1;
        
    end
    
    
    % Remove duplicated triangles stuck together
    % Attention à ne pas créer de trous dans le maillage !
    T_sort = sort(T_cur,2);
    [~,j] = unique(T_sort,'rows');
    g = setdiff(1:size(T_sort,1),j);
    f = zeros(0,1);
    
    for p = 1:numel(g)
       
        f = cat(1,f,find(ismember(T_sort,T_sort(g(p),:),'rows')));
        
    end
    
    T_cur(f,:) = [];
    N_cur(f,:) = [];    
    
end


end % multiresolution_mesh9


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