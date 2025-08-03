function [] = convex_hull_animation()
%% Jarvis / gift wrapping algorithm
%
% Author : nicolas.douillet9 (at) gmail.com, 2020-2025.


load('../../data/Random_unit_ball_128_pts.mat');
on_vtx_id = false(size(V,1),1);


%%  ------------------- Display parameters --------------------- %
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

title_text = {'Convex hull construction of a random spherical point set', 'Via gift wrapping / Jarvis algorithm'};
filename = 'convex_hull_animation.gif';


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

% Initial point set only (red / yellow)
plot3(V(~on_vtx_id,1),V(~on_vtx_id,2),V(~on_vtx_id,3),strcat(vtx_color_in,vertex_marker),'MarkerSize',vertex_size,'MarkerEdgeColor',vtx_color_in,'MarkerFaceColor',vtx_color_in,'LineWidth',3);
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
nb_vtx = size(V,1); % nb vertices

% 1st vertex : the one minimum z
[~,v1] = mink(V(:,3),1);

% 2nd vertex : the one which build an edge with minimum (positive) dot product with k [0 0 1] vector
d = dot(V-V(v1,:),repmat([0 0 1],[size(V,1),1]),2);
[~,v2] = mink(d,2);   
v2 = setdiff(v2,v1); % remove dot prod with null vector

curr_edge = [v1 v2]; % initial edge
on_vtx_id([v1,v2],1) = true;
edg_list = curr_edge;
proc_edg_list = zeros(0,2); % processed edge list
T = zeros(0,3);
curr_edge_id = 1;


while ~isempty(curr_edge)
    
    [new_tgl,nxt_vtx_id] = find_nxt_vertex(curr_edge,V,nb_vtx,epsilon); % next vertices
    proc_edg_list = cat(1,proc_edg_list,curr_edge);
    curr_edge_id = curr_edge_id + 1;
    
    % Non already existing new triangles
    idx2keep = ~ismember(sort(new_tgl,2),sort(T,2),'rows');
    new_tgl = new_tgl(idx2keep,:);
    nxt_vtx_id = nxt_vtx_id(idx2keep);    
    on_vtx_id(nxt_vtx_id) = true;
    
    if ~isempty(new_tgl)
        
        T = cat(1,T,new_tgl); % add triangles        
        
        
        plot3(V(~on_vtx_id,1),V(~on_vtx_id,2),V(~on_vtx_id,3),strcat(vtx_color_in,vertex_marker),'MarkerSize',vertex_size,'MarkerEdgeColor',vtx_color_in,'MarkerFaceColor',vtx_color_in,'LineWidth',edge_width), hold on;                 
        plot3(V(on_vtx_id,1), V(on_vtx_id,2), V(on_vtx_id,3), strcat(vtx_color_on,vertex_marker),'MarkerSize',vertex_size,'MarkerEdgeColor',vtx_color_on,'MarkerFaceColor',vtx_color_on,'LineWidth',edge_width);                
        t = trisurf(T,V(:,1),V(:,2),V(:,3),'LineWidth',2); shading flat;
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
        
        % Create and add the two new edges
        new_edg1 = cat(2,repmat(curr_edge(1),[numel(nxt_vtx_id),1]),nxt_vtx_id');
        new_edg2 = cat(2,repmat(curr_edge(2),[numel(nxt_vtx_id),1]),nxt_vtx_id');                
        edg_list = cat(1,edg_list,new_edg1,new_edg2);
        
    end
        
    % Compute new edge index to process
    if curr_edge_id < 1 + size(edg_list,1)
        
        nxt_edg = edg_list(curr_edge_id,:);
        
        if ~ismember(sort(nxt_edg,2),sort(proc_edg_list,2),'rows')
            
            curr_edge = nxt_edg;
            
        end
        
    else
        
        curr_edge = [];
        
    end
    
end


plot3(V(~on_vtx_id,1),V(~on_vtx_id,2),V(~on_vtx_id,3),strcat(vtx_color_in,vertex_marker),'MarkerSize',vertex_size,'MarkerEdgeColor',vtx_color_in,'MarkerFaceColor',vtx_color_in,'LineWidth',edge_width), hold on;
plot3(V(on_vtx_id,1), V(on_vtx_id,2), V(on_vtx_id,3), strcat(vtx_color_on,vertex_marker),'MarkerSize',vertex_size,'MarkerEdgeColor',vtx_color_on,'MarkerFaceColor',vtx_color_on,'LineWidth',edge_width), hold on;
t = trisurf(T,V(:,1),V(:,2),V(:,3),'LineWidth',2); shading faceted, hold on;
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


end % convex_hull_animation


%% find_nxt_vertex subfunction
function [new_tgl,nxt_vtx_id] = find_nxt_vertex(edge, V, nb_vtx, epsilon)
% find_nxt_vertex : function to find the next vertices in the algorithm
% given a current edge. These candidate new vertices are the ones which
% build a triangle which has all the other vertices of the point set on one
% same unique side.


new_tgl = zeros(0,3);
nxt_vtx_id = [];

% Candidate triangles
vtx_id2test = setdiff(1:nb_vtx,edge);
tglist = cat(2,repmat(edge,[nb_vtx-2,1]),vtx_id2test');
G = cell2mat(cellfun(@(r) mean(V(r,:),1),num2cell(tglist,2),'un',0));

% Normal vectors
tgl_normals = cross(V(tglist(:,3),:)-V(tglist(:,1),:),V(tglist(:,2),:)-V(tglist(:,1),:),2);


for i = 1:numel(vtx_id2test)
    
    cnd_vtx_vect = V(vtx_id2test,:) - repmat(G(i,:),[nb_vtx-2,1]); 
    cnd_vtx_vect = cnd_vtx_vect ./ sqrt(sum(cnd_vtx_vect.^2,2));
    
    dot_prod = dot(cnd_vtx_vect,tgl_normals,2);
    dot_prod(abs(dot_prod) < epsilon) = 0;
    sgn_dot_prod = sign(dot_prod);  
    
    if isequal(sgn_dot_prod,abs(sgn_dot_prod)) % native outward oriented triangles concatenation
        
        new_tgl = cat(1,new_tgl,tglist(i,:));
        nxt_vtx_id = cat(2,nxt_vtx_id,vtx_id2test(i));
        
    elseif isequal(sgn_dot_prod,-abs(sgn_dot_prod)) % native inward oriented triangles reorientation + concatenation
        
        new_tgl = cat(1,new_tgl,fliplr(tglist(i,:)));
        nxt_vtx_id = cat(2,nxt_vtx_id,vtx_id2test(i));
        
    end
    
end


end % find_new_vertex