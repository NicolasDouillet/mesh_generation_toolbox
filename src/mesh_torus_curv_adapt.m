function [V, T] = mesh_torus_curv_adapt(R, r)
% mesh_torus_curv_adapt : function to create mesh adaptated to the curvature for the torus.
%
% Author & support : nicolas.douillet (at) free.fr, 2023.
%
%
%%% Input arguments :
%
% - R > r > 0 : positive real scalar double, the torus large radius.
% - r > 0 :     positive real scalar double, the torus small radius.
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


% Input parsing
assert(r > 0, 'Torus small radius must be a real positive number');
assert(R > 0, 'Torus great radius must be a real positive number');
assert(R > r, 'Torus big radius must be greater than torus small radius.');

%                                  ____________________________________________
%                                 |                                            |
nb_samples = 1 + floor(R/r); % <= | modify the sampling rate HERE if necessary |
%                                 |____________________________________________|

% Build icosahedron
phi_n = 0.5*(1+sqrt(5));
centre_angle = 2*asin(1/sqrt(phi_n*sqrt(5)));

Rmz = @(theta)[cos(theta) -sin(theta) 0;...
               sin(theta)  cos(theta) 0;...
               0          0           1];

% 1st equilateral triangle
V1 = [0 0 1]';
V2 = [sin(centre_angle) 0 cos(centre_angle)]';
V3 = Rmz(0.4*pi)*V2;

% Lower base triangle with /O symetry
V12 = -V1;
V11 = -V2;
V10 = -V3;

% (12) Vertices set coordinates vector
V4 = Rmz(0.4*pi)*V3;
V5 = Rmz(0.8*pi)*V3;
V6 = Rmz(1.2*pi)*V3;
V9 = Rmz(0.4*pi)*V10;
V8 = Rmz(0.8*pi)*V10;
V7 = Rmz(1.2*pi)*V10;

V = [V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12]';

T = [3 2 8; % centre belt (10 triangles)
     3 8 7;
     4 3 7;
     4 7 11;
     5 4 11;
     5 11 10;
     6 5 10;
     6 10 9;
     2 6 9;
     2 9 8];

% Remove icosahedron top and bottom vertices
V = V(2:end-1,:);
T = T - 1;
     

% Sample triangles
for k = 1:size(T,1)
       
    [Vt,Tt] = mesh_triangle(V(T(k,1),:)',V(T(k,2),:)',V(T(k,3),:)',nb_samples);
            
    T = cat(1,T,Tt+size(V,1));
    V = cat(1,V,Vt);    

end

% Remove original vertex and triangle set from the icosahedron
if nb_samples > 1
    
    V = V(11:end,:);
    T = T(11:end,:);
    T = T - 10;
    
end


% Symetric lower part
V(:,3) = V(:,3) - min(V(:,3));
T = cat(1,T,fliplr(T+size(V,1)));
V = cat(1,V,cat(2,V(:,1),V(:,2),-V(:,3)));


% Stretch the box here in X/Y/Z
V(:,3) = V(:,3)/max(V(:,3));
V(:,3) = pi*(V(:,3)/12-min(V(:,3)));
V(:,1:2) = r*V(:,1:2) ./ sqrt(sum(V(:,1:2).^2,2));


Rmy = @(theta)[cos(theta) 0 -sin(theta);
               0          1  0;
               sin(theta) 0  cos(theta)];

V(:,1) = R + V(:,1) - r;
Z = V(:,3);


% Crescent creation
for k = 1:size(V,1)
    
    V(k,:) = (Rmy(V(k,3))*cat(2,V(k,1:2),0)')';
    
end

V(:,3) = V(:,3) + r*sin(Z);
V(:,1) = V(:,1) + r*cos(Z);
    

% Concatenations
T = cat(1,T,T+size(V,1));
V = cat(1,V,(Rmy(pi/6)*V')');
T = cat(1,T,T+size(V,1));
V = cat(1,V,(Rmy(pi/3)*V')');
T = cat(1,T,T+size(V,1),T+2*size(V,1));
V = cat(1,V,(Rmy(2*pi/3)*V')',(Rmy(4*pi/3)*V')');

% Switch coordinates to have (Oz) as the main symetry axis
V = cat(2,V(:,1),V(:,3),V(:,2));


% Remove duplicated vertices
[V,~,n] = uniquetol(V,1e3*eps,'ByRows',true);
T = n(T);

% Remove duplicated triangles
T_sort = sort(T,2);
[~,idx,~] = unique(T_sort,'rows','stable');
T = T(idx,:);


end % mesh_torus_curv_adapt