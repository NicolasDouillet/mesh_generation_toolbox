function [V, T] = mesh_ovoid(nb_samples)
%% mesh_ovoid : function to mesh an ovoid following
% Hügelschäffer egg curve parameteric equation.
%
% Author : nicolas.douillet (at) free.fr, 2021-2024.
%
%
% Input arguments :
%
% - nb_samples > 1 : positive integer scalar double, the number of samples / levels 
%                    along one side of the ovoïd curve, from top to bottom.
%
%
% Output arguments :
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
% About / other informations
%
% Ovoid longitudinal axis is (Oz).
% Triangles / normals are not coherently oriented.


%% Body
angle_step = pi/nb_samples;

% Z axis rotation matrix
Rmz = @(theta)[cos(theta) -sin(theta) 0;
               sin(theta)  cos(theta) 0;
               0           0          1];
      
angle_vect = angle_step:angle_step:2*pi;     

% Dafault parameter values
a = 6;
b = 4;
d = 1;

% Hügelschäffer egg curve parameteric equation
Z = (sqrt(a^2 - d^2*sin(angle_vect).^2) + d*cos(angle_vect)).*cos(angle_vect);
X = b*sin(angle_vect);
Y = zeros(1,numel(X));

U = cat(1,X,Y,Z);
V = U';          

for theta = angle_vect
    
    R = (Rmz(theta)*U)';
    V = cat(1,V,R);
    
end


S1 = numel(X);       % number of colums
S2 = 2*nb_samples+1; % number of lines

% Theoritical number of quadrangles : S1*S2 
% Theoritical number of triangles : 2*S1*S2,
% Mais aux pôles les triangles sont confondus
% -> bien s'assurer qu'ils sont tous orientés dans le même sens
T = build_triangulation(S1,S2,nb_samples);      
           
% Remove duplicated vertices
tol = 1e3*eps;
[V,T] = remove_duplicated_vertices(V,T,tol);

% Remove duplicated triangles
T = remove_duplicated_triangles(T);


end % mesh_ovoid


%% build_triangulation subfunction
function T = build_triangulation(S1, S2, nb_samples)


r1 = cat(2,1,repelem(2:S1-1,2),S1);
r1 = reshape(r1,[2,S1-1])';
R1 = cat(2,r1,(1:S1-1)'+S1); % 1st triangle row indices
R1 = cat(1,R1,[2*nb_samples 2*nb_samples+1 1]); % add last triangle to connect
R1(1+floor(0.5*size(R1,1)):end,:) = fliplr(R1(1+floor(0.5*size(R1,1)):end,:));

r2 = cat(2,1+S1,repelem(2+S1:2*S1-1,2),2*S1);
r2 = reshape(r2,[2,S1-1])';
R2 = cat(2,(2:S1)',fliplr(r2)); % 2nd triangle row indices
% R2(1+floor(0.5*size(R2,1)):end,:) = fliplr(R2(1+floor(0.5*size(R2,1)):end,:));

T = repmat(cat(1,R1,R2),[S2-2,1]);
T = T + floor(0.5*S1)*repelem((0:S2-3)',2*(S1-1)+1,3);


end % build_triangulation