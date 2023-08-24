function [V, T] = mesh_torus_specific(R, r)
% mesh_torus_specific : function to mesh a torus while allowing to specifically choose
% the number of samples in the two directions, longitude and latitude.
%
% Author & support : nicolas.douillet (at) free.fr, 2023.
%
%
%%% Input arguments :
%
% - R > r > 0 : positive integer scalar double, the torus large radius.
% - r > 0     : positive ODD scalar double, the torus small radius.
%
%
%%% Output arguments :
%
%        [| | |]
% - V_ = [X Y Z], real matrix double, the output point set, size(V) = [nb_vertices,3]
%        [| | |]
%                 with nb_vertices is the same number as the number of vertices in P.
%
%       [|  |  |]
% - T = [i1 i2 i3], positive integer matrix double, the triangulation, size(T) = [nb_triangles,3].
%       [|  |  |]
%
%
%%% About / other informations
%
% Torus is centered on the origin, [0 0 0], and has (0z) for axis.
% Triangles / normals are coherently oriented and facing outward.


a = 1; % mesh step
M = regular_2D_trigrid(R,r,a);


% Indexation triangles
T = index_mesh_triangles(M);
V = cat(2,reshape(M,[numel(M(:,:,1)),2]),zeros(numel(M(:,:,1)),1));

V(:,1) = V(:,1)/max(V(:,1));
V(:,1) = V(:,1)*2*pi*(1+0.5*a/(R-1)); 
V(:,2) = V(:,2)/max(V(:,2));
V(:,2) = V(:,2)*2*pi;

X = (R + r*sin(V(:,2))).*cos(V(:,1));
Y = (R + r*sin(V(:,2))).*sin(V(:,1));
Z = r*cos(V(:,2));
V = cat(2,X,Y,Z);


% Remove duplicated vertices
[V,~,n] = uniquetol(V,1e3*eps,'ByRows',true);
T = n(T);

% Remove duplicated triangles
T_sort = sort(T,2);
[~,idx,~] = unique(T_sort,'rows','stable');
T = T(idx,:);


end % mesh_torus_specific


% index_mesh_triangles subfunction
function T = index_mesh_triangles(M)
%
% Author & support : nicolas.douillet (at) free.fr, 2023.


H = size(M,2);
W = size(M,1);

T = [];


for i = 1:H
   
    for j = 1:W-1
        
        if mod(i,2) == 1
           
            if i < H
                
                T = cat(1,T,[(i-1)*W+j,(i-1)*W+j+1,i*W+j]);
                
            end

            if i > 1
                
                T = cat(1,T,[(i-1)*W+j,(i-1)*W+j+1,(i-2)*W+j]);
                
            end
                                    
        else % if mod(i,2) == 0
            
            
            if i < H
                
                T = cat(1,T,[(i-1)*W+j,(i-1)*W+j+1,i*W+j+1]);
            
            end
            
            T = cat(1,T,[(i-1)*W+j,(i-1)*W+j+1,(i-2)*W+j+1]);
                
        end
        
    end
    
end


end % index_mesh_triangles


% regular_2D_trigrid subfunction
function M = regular_2D_trigrid(W, H, a)
%
% Author & support : nicolas.douillet (at) free.fr, 2023.


h = 0.5*a*sqrt(3); % equilateral triangle height

% X Y coordinates computation and concatenation
Vx = repmat(linspace(0,W-1,W),[H,1]); 
Vx(mod(1:H,2)==0,:) = Vx(mod(1:H,2)==0,:) + 0.5*a;
Vy = repmat(h*linspace(0,H-1,H)',[1,W]);

M = cat(3,Vx',Vy');


end % regular_2D_trigrid