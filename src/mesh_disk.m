function [V, T] = mesh_disk(r, edg_nb_smpl)
% mesh_disk : function to mesh a disk.
%
% Author & support : nicolas.douillet (at) free.fr, 2023.
%
%
%%% Input arguments :
%
% - r : positive real scalar double, the disk radius.
%
% - edg_nb_smpl : positive integer scalar double, one sixth of the number of samples on the disk perimeter.
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
% Disk is centered on the origin, [0 0 0], and has (Oz) for normal axis.
% Triangles / normals are coherently oriented.


M1 = [0.5*sqrt(3)  0.5 0]';
M2 = [0.5*sqrt(3) -0.5 0]';

[V1,T1] = mesh_triangle(zeros(3,1),M1,M2,edg_nb_smpl);


% Disk sector transformation
% Intersection points with the vertical line
X = zeros(size(V1));
X(:,1) = 0.5*sqrt(3);
X(:,3) = zeros(size(V1,1),1);
X(:,2) = V1(:,2) + (0.5*sqrt(3) - V1(:,1)).*V1(:,2)./V1(:,1);

N1 = sqrt(sum(V1.^2,2));
Nx = sqrt(sum(X.^2,2));
norm_pct = N1./Nx;

% Update norm = distance ratio OM/OX
V1 = r*V1.*norm_pct./N1;
V1(1,:) = zeros(1,3);

% Duplication by rotation
Rmz = @(theta)[cos(theta) -sin(theta) 0;
               sin(theta)  cos(theta) 0;
               0           0          1];

V2 = (Rmz(pi/3)*V1')';
V3 = (Rmz(2*pi/3)*V1')';
V4 = cat(1,V1,V2,V3);
V5 = (Rmz(pi)*V4')';

V = cat(1,V4,V5);
T = cat(1,T1,T1+size(V1,1),T1+2*size(V1,1),T1+3*size(V1,1),T1+4*size(V1,1),T1+5*size(V1,1));

% Duplicata removal
[V,~,n] = uniquetol(V,'ByRows',true);
T = n(T);


end % mesh_disk