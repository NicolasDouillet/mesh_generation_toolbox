function [] = quick_hull_animation()
%% Simplex / divide and conquer algorithm
%
% Author : nicolas.douillet9 (at) gmail.com, 2020-2025.


load('../../data/Random_unit_ball_128_pts.mat');
V_in = V;


%% -------------------- Display parameters --------------------- %
time_lapse = 0.1;      % animation time lapse; default value : 0.5
bckgrd_clr = [0 0 0];  % background color; default : [0 0 0] (black)
text_color = [1 1 1];  % title color; default : [1 1 1] (white)
vtx_color_in  = 'r';   % Inside vertex color; default : 'r' (red)
vtx_color_on  = 'g';   % On convex hull vertex color; default : 'g' (green)
face_color = [0 1 1];  % face color; default :   [0 1 1] (cyan)
edge_color = [0 0 0];  % edge color; default :   [0 1 1] (cyan)
edge_width = 2;        % edge width; default : 2
vertex_marker = '+';   % vertex marker; default : '+' spare : 's'
vertex_size = 8;       % vertex size; default : 6
el = 5;                % view elevation

title_text = {'Animation of the convex hull construction steps of a random 3D point set', 'Computed with quick hull / divide and conquer algorithm'};
filename = 'quick_hull.gif';


% -------- Display settings -------- %
h = figure;
set(h,'Position',get(0,'ScreenSize'));
set(gcf,'Color',bckgrd_clr);

drawnow;
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif', 'Loopcount',Inf,'DelayTime',time_lapse);
clf;

plot3(V_in(:,1),V_in(:,2),V_in(:,3),strcat(vtx_color_in,vertex_marker),'MarkerSize',vertex_size,'MarkerEdgeColor',vtx_color_in,'MarkerFaceColor',vtx_color_in,'LineWidth',3), hold on;
set(gca,'Color',bckgrd_clr);
axis equal tight;
axis off;
view(0,el);
% title(title_text,'Color',text_color,'FontSize',16);

drawnow;
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',2*time_lapse);
clf;


%% Body
coeff = 1;
epsilon = coeff*eps; % floating point tolerance error
nb_vtx = size(V_in,1);

f_Xmin = find(V_in(:,1) == min(V_in(:,1)));
f_Xmin = f_Xmin(1,1);
f_Ymin = find(V_in(:,2) == min(V_in(:,2)));
f_Ymin = f_Ymin(1,1);
f_Zmin = find(V_in(:,3) == min(V_in(:,3)));
f_Zmin = f_Zmin(1,1);
f_Xmax = find(V_in(:,1) == max(V_in(:,1)));
f_Xmax = f_Xmax(1,1);
f_Ymax = find(V_in(:,2) == max(V_in(:,2)));
f_Ymax = f_Ymax(1,1);
f_Zmax = find(V_in(:,3) == max(V_in(:,3)));
f_Zmax = f_Zmax(1,1);

bounding_box = [f_Xmin f_Xmax f_Ymin f_Ymax f_Zmin f_Zmax];
hull_vtx_id = unique(bounding_box);

if numel(hull_vtx_id) >= 4
   
    hull_vtx_id = hull_vtx_id(1,1:4);
    T = combnk(1:4,3);
        
elseif numel(hull_vtx_id) == 3
    
    T = 1:3; 
    
else % if numel(hull_vtx_id) < 3
    
    error('Colinear input point set.');
    
end

N = face_normals(V_in(hull_vtx_id,:),T,'norm');
    
% Avoid initial flat tetrahedron cases (octahedron example)
while isequal(N,repmat(N(1,:),[size(N,1),1]))
    
    hull_vtx_id = [];
    
    while numel(hull_vtx_id) < 4
        
        bounding_box = randi(nb_vtx,1,4);
        hull_vtx_id = unique(bounding_box);        
        
    end
        
    N = face_normals(V_in(hull_vtx_id,:),T,'norm');
    
end

nb_t = 4;
hull_vtx_id = hull_vtx_id(1:nb_t);


T = nchoosek(hull_vtx_id,3);
G = mean(V_in(hull_vtx_id,:),1);
N = face_normals(V_in,T,'norm');

% Orient normals outward
Gt = cell2mat(cellfun(@(r) mean(V_in(r,:),1),num2cell(T,2),'un',0));
orientation = sign(dot(N,Gt-repmat(G,[nb_t,1]),2));

if ~isequal(orientation,ones(nb_t,1)) && ~isequal(orientation,zeros(nb_t,1))
    N = N.*orientation;
    T(orientation < 0,:) = fliplr(T(orientation < 0,:));
end

[V_out,T] = remove_inside_pts(V_in,T,epsilon);
on_vtx_id = false(size(V_out,1),1);
on_vtx_id(unique(T(:)),1) = true;


Vx = setdiff(V_in,V_out(on_vtx_id,:),'rows');

plot3(Vx(:,1),Vx(:,2),Vx(:,3),strcat(vtx_color_in,vertex_marker),'MarkerSize',vertex_size,'MarkerEdgeColor',vtx_color_in,'MarkerFaceColor',vtx_color_in,'LineWidth',3), hold on;
plot3(V_out(on_vtx_id,1), V_out(on_vtx_id,2), V_out(on_vtx_id,3), strcat(vtx_color_on,vertex_marker),'MarkerSize',vertex_size,'MarkerEdgeColor',vtx_color_on,'MarkerFaceColor',vtx_color_on,'LineWidth',edge_width);   
t = trisurf(T,V_out(:,1),V_out(:,2),V_out(:,3),'LineWidth',2); shading flat;
set(t,'FaceColor',face_color,'EdgeColor',edge_color);
set(gca,'Color',bckgrd_clr);
axis equal tight;
axis off;
alpha(0.5);

view(0,el);
% title(title_text,'Color',text_color,'FontSize',16);

drawnow;
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',time_lapse);
clf;


nb_new_tgl = true;

while nb_new_tgl
    
    curr_tgl_id = 1;
    nb_new_tgl = 0;
    
    while curr_tgl_id < 1 + size(T,1)
                
        [T,N,new_vtx_id] = grow_tetrahedron(V_out,T,N,curr_tgl_id,epsilon);
        on_vtx_id(new_vtx_id,1) = true;
        
        Vx = setdiff(V_in,V_out(on_vtx_id,:),'rows');
       
        plot3(Vx(:,1),Vx(:,2),Vx(:,3),strcat(vtx_color_in,vertex_marker),'MarkerSize',vertex_size,'MarkerEdgeColor',vtx_color_in,'MarkerFaceColor',vtx_color_in,'LineWidth',3), hold on;
        plot3(V_out(on_vtx_id,1), V_out(on_vtx_id,2), V_out(on_vtx_id,3), strcat(vtx_color_on,vertex_marker),'MarkerSize',vertex_size,'MarkerEdgeColor',vtx_color_on,'MarkerFaceColor',vtx_color_on,'LineWidth',edge_width);   
        t = trisurf(T,V_out(:,1),V_out(:,2),V_out(:,3),'LineWidth',2); shading flat;
        set(t,'FaceColor',face_color,'EdgeColor',edge_color);
        set(gca,'Color',bckgrd_clr);
        axis equal tight;
        axis off;
        alpha(0.5);
        
        view(0,el);
        % title(title_text,'Color',text_color,'FontSize',16);
        
        drawnow;
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',time_lapse);
        clf;
        
        
        if new_vtx_id % effective grow with new triangles
                              
            nb_new_tgl = nb_new_tgl + 2;
            edg_list = query_edg_list(T,'sorted');
            i = 1;
            
            while i < 1 + size(edg_list,1)
                
                tgl_pair_id = cell2mat(find_triangle_indices_from_edg_list(T,edg_list(i,:)));                
                isconcave = detect_concavity(V_out,T,N,tgl_pair_id,epsilon);
                
                if isconcave
                    
                    [T,N,edg_list] = flip_two_ngb_triangles(tgl_pair_id,T,V_out,N,edg_list);
                    
                    Vx = setdiff(V_in,V_out(on_vtx_id,:),'rows');
                    
                    plot3(Vx(:,1),Vx(:,2),Vx(:,3),strcat(vtx_color_in,vertex_marker),'MarkerSize',vertex_size,'MarkerEdgeColor',vtx_color_in,'MarkerFaceColor',vtx_color_in,'LineWidth',3), hold on;
                    plot3(V_out(on_vtx_id,1), V_out(on_vtx_id,2), V_out(on_vtx_id,3), strcat(vtx_color_on,vertex_marker),'MarkerSize',vertex_size,'MarkerEdgeColor',vtx_color_on,'MarkerFaceColor',vtx_color_on,'LineWidth',edge_width);   
                    t = trisurf(T,V_out(:,1),V_out(:,2),V_out(:,3),'LineWidth',2); shading flat;
                    set(t,'FaceColor',face_color,'EdgeColor',edge_color);
                    set(gca,'Color',bckgrd_clr);
                    axis equal tight;
                    axis off;
                    alpha(0.5);
                    
                    view(0,el);
                    % title(title_text,'Color',text_color,'FontSize',16);
                    
                    drawnow;
                    frame = getframe(h);
                    im = frame2im(frame);
                    [imind,cm] = rgb2ind(im,256);
                    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',time_lapse);
                    clf;
                    
                    
                else
                    
                    i = i + 1;
                    
                end
                
            end
            
            [V_out,T] = remove_inside_pts(V_out,T,epsilon);            
            on_vtx_id = false(size(V_out,1),1);
            on_vtx_id(unique(T(:)),1) = true;
            curr_tgl_id = curr_tgl_id - 1;
              
        else
            
            nb_new_tgl = 0;
            
        end
        
        curr_tgl_id = curr_tgl_id + 1;
        
    end       
    
end

Vx = setdiff(V_in,V_out(on_vtx_id,:),'rows');

plot3(Vx(:,1),Vx(:,2),Vx(:,3),strcat(vtx_color_in,vertex_marker),'MarkerSize',vertex_size,'MarkerEdgeColor',vtx_color_in,'MarkerFaceColor',vtx_color_in,'LineWidth',3), hold on;
plot3(V_out(:,1),V_out(:,2),V_out(:,3),strcat(vtx_color_on,vertex_marker),'MarkerSize',vertex_size,'MarkerEdgeColor',vtx_color_on,'MarkerFaceColor',vtx_color_on,'LineWidth',3), hold on;
t = trisurf(T,V_out(:,1),V_out(:,2),V_out(:,3),'LineWidth',2); shading flat;
set(t,'FaceColor',face_color,'EdgeColor',edge_color);
set(gca,'Color',bckgrd_clr);
axis equal tight;
axis off;
alpha(0.5);


% Roll 360° view point
angle_step = 5;

for phi = angle_step:angle_step:360-angle_step
    
    view(phi,el);
    % title(title_text,'Color',text_color,'FontSize',16);
    
    drawnow;
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',time_lapse);
    
end

close(h);


end % quick_hull


%% grow_tetrahedron subfunction
function [T, N, new_vtx_id] = grow_tetrahedron(V, T, N, tgl_id, epsilon)
% grow_tetrahedron : function to find one
% new vertex belonging to the convex hull
% for one given triangle, to create the three
% newborn triangles and erase the previous one.


new_vtx_id = [];
nb_vtx = size(V,1);

d = dot(repmat(N(tgl_id,:),[nb_vtx,1]),V,2);

if find(abs(d) > epsilon)
    
    f = find(d == max(d));    
    
    if f & ~ismember(f,unique(T(:))) % prevent from creating non manifold stuffs
        
        f = f(1,1);
        
        new_tgl1 = cat(2,T(tgl_id,1:2),f);
        new_tgl2 = cat(2,T(tgl_id,2:3),f);
        new_tgl3 = cat(2,T(tgl_id,3),T(tgl_id,1),f);
                
        % Add 3 new triangles and face normals
        T = cat(1,T,new_tgl1,new_tgl2,new_tgl3);
        new_face_normals = face_normals(V,T(end-2:end,:),'norm');
        N = cat(1,N,new_face_normals);           
        
        % Remove one triangle and its normal
        T(tgl_id,:) = [];        
        N(tgl_id,:) = [];
        new_vtx_id = f;
        
    end
    
end


end % grow_tetrahedron


%% detect_concavity subfunction
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


%% flip_two_ngb_triangles subfunction
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


%% remove_inside_pts subfunction
function [V_out, T] = remove_inside_pts(V_in, T, epsilon)
% remove_inside_pts : function to remove points inside the convex hull
% during its computational iterations to save cpu time.


in_vtx_set_id = find(isin3Dconvexset(V_in,T,V_in,epsilon));

if ~isempty(in_vtx_set_id)
    
    [V_out,T] = remove_vertices(in_vtx_set_id,V_in,T,'index');
    
else
    
    V_out = V_in;
    
end


end % remove_inside_pts