function [V, T] = mesh_reuleaux_tetrahedron(sample_step, option_display)
% mesh_reuleaux_tetrahedron : function to compute, display, and save a meshed Reuleaux tetrahedron. 
%
% Author & support : nicolas.douillet (at) free.fr, 2017-2023.
%
%
% Syntax
%
% mesh_reuleaux_tetrahedron;
% mesh_reuleaux_tetrahedron(sample_step);
% mesh_reuleaux_tetrahedron(sample_step, option_display);
% [V, T] = mesh_reuleaux_tetrahedron(sample_step, option_display);
%
%
% Description
%
% mesh_reuleaux_tetrahedron computes and displays the meshed Reuleaux
% tetrahedron included in the unit sphere, and which each
% edge is sampled in 8.
%
% mesh_reuleaux_tetrahedron(sample_step) uses sample_step steps.
%
% mesh_reuleaux_tetrahedron(sample_step, option_display)
% displays the result when option_display is logical *true/1, and doesn't when it is
% logical false/0.
%
% [V, T] = mesh_reuleaux_tetrahedron(sample_step, option_display) stores the resulting
% vertices coordinates in the array V, and the corresponding triplet indices list in the array T.
% 
%
% See also
%
% <https://fr.mathworks.com/help/matlab/ref/mesh.html?s_tid=srchtitle mesh> | 
% <https://fr.mathworks.com/help/matlab/ref/trimesh.html?searchHighlight=trimesh&s_tid=doc_srchtitle trimesh>
%
%
% Input arguments
%
% - sample_step : positive integer scalar, power of 2.
%
% - option_display : logical *true (1) / false (0).
%
%
% Output arguments
%
%     [ |  |  |]
% - V [Vx Vy Vz] : real matrix double, the point set. Size = [nb_vertices,3].
%     [ |  |  |]
%
%     [ |  |  |]
% - T [T1 T2 T3] : positive integer matrix double, the triangulation. Size = [nb_triangles,3].
%     [ |  |  |]


% TODO : nb_steps rotations du bord �chantillon� autour de l'axe d�fini par le sommet sup�rieur et le sommet oppos�, puis triangulation. 


% Input parsing
assert(nargin < 3,'Too many input arguments.');

if nargin > 0
    assert(isnumeric(sample_step) && sample_step == floor(sample_step) && sample_step > 0,'sample_step parameter must be a positive integer.');    
    if nargin  > 1                                
            assert(islogical(option_display) || isnumeric(option_display),'option_display parameter type must be either logical or numeric.');                    
    else                
        option_display = true;        
    end    
else
    sample_step = 8;    
    option_display = true;    
end


% Body
%  Summits of original tetrahedron, living in the sphere S(O,1)
V1 = [0 0 1];
V2 = [2*sqrt(2)/3 0 -1/3];
V3 = [-sqrt(2)/3 sqrt(6)/3 -1/3];
V4 = [-sqrt(2)/3 -sqrt(6)/3 -1/3];

edge_length = norm(V1-V2); %  = 2*sqrt(6)/3

[V123, T] = mesh_triangle(V2', V1', V3', sample_step);
V_flat = V123; 

D123 = sqrt(sum((V123 - V4).^2,2)); % distance matrix
V123 = edge_length*(V123 - V4) ./ repmat(D123, [1 3]) + repmat(V4, [size(V123,1), 1]); % "inflated triangle" / opposite vertex V4

be = 1:sample_step+1;                       % bottom edge index vector
re = cumsum(sample_step+1:-1:0);            % right edge index vector
le = cat(2,1,1+cumsum(sample_step+1:-1:2)); % left edge index vector

M34 = [-sqrt(2)/3 0 -1/3]; % middle of [V3;V4] segment
D34 = sqrt(sum((V_flat(le,:) - M34).^2,2));
V_flat(le,:) = sqrt(2)*(V_flat(le,:) - M34) ./ repmat(D34, [1 3]) + repmat(M34, [size(V_flat(le,:),1), 1]);
V123(le,:) = V_flat(le,:);

M24 = [1/3/sqrt(2) -1/sqrt(6) -1/3]; % middle of [V2;V4] segment
D24 = sqrt(sum((V_flat(re,:) - M24).^2,2));
V_flat(re,:) = sqrt(2)*(V_flat(re,:) - M24) ./ repmat(D24, [1 3]) + repmat(M24, [size(V_flat(re,:),1), 1]);
V123(re,:) = V_flat(re,:);

M14 = [-1/3/sqrt(2) -1/sqrt(6) 1/3]; % middle of [V1;V4] segment
D14 = sqrt(sum((V_flat(be,:) - M14).^2,2));
V_flat(be,:) = sqrt(2)*(V_flat(be,:) - M14) ./ repmat(D14, [1 3]) + repmat(M14, [size(V_flat(be,:),1), 1]);
V123(be,:) = V_flat(be,:);


% Tetrahedron faces rotations
Rmy = @(theta) [cos(theta) 0 -sin(theta);
                0          1  0;
                sin(theta) 0  cos(theta)];

Rmz = @(theta) [cos(theta) -sin(theta) 0;
                sin(theta)  cos(theta) 0;
                0           0          1];
                        
V134 = (Rmz(2*pi/3)*V123')';
V142 = (Rmz(2*pi/3)*V134')';
V234 = (Rmy(-acos(-1/3))*Rmz(pi/3)*V142')';
            
V = [V123; V134; V142; V234];

% Triplet indices list 
T = [T;
     T+repmat(size(V123,1),   [size(T,1) size(T,2)]);...
     T+2*repmat(size(V123,1), [size(T,1) size(T,2)]);...
     T+3*repmat(size(V123,1), [size(T,1) size(T,2)])];
    

% Duplicated vertices removal (from the edges)
[V,T] = remove_duplicated_vertices(V,T);


TRI = triangulation(T, V(:,1), V(:,2), V(:,3));

%  Display
if option_display
    
    figure;    
    trimesh(TRI), hold on;    
    colormap([0 0 1]);    
    axis square, axis equal, axis tight;
    
end


end % mesh_reuleaux_tetrahedron


% remove_duplicated_vertices subfunction
function [V_out, T_out] = remove_duplicated_vertices(V_in, T_in)

tol = 1e4*eps;
[V_out,~,n] = uniquetol(V_in,tol,'ByRows',true);
T_out = n(T_in);

end % remove_duplicated_vertices