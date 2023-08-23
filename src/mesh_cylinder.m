function [V, T] = mesh_cylinder(r, h, edg_nb_smpl)
%
% Author & support : nicolas.douillet (at) free.fr, 2023.
%
%
% - r : the cylinder radius.
% - h : the cylinder height.
% - edg_nb_smpl : one sixth of the number of samples on the disk perimeter.
%
% Cylinder is centred on the origin, [0 0 0], and has (Oz) for axis.


M1 = [0.5*r*sqrt(3)  0.5*r  0.5*h];
M2 = [0.5*r*sqrt(3) -0.5*r  0.5*h];
M3 = [0.5*r*sqrt(3) -0.5*r -0.5*h];
M4 = [0.5*r*sqrt(3)  0.5*r -0.5*h];

[V1,T1] = mesh_quadrangle(M1, M2, M3, M4, edg_nb_smpl);

% Cylinder sector transformation
% Intersection points with the vertical plane
X = zeros(size(V1));
X(:,1) = 0.5*r*sqrt(3);
X(:,2) = V1(:,2) + (0.5*r*sqrt(3) - V1(:,1)).*V1(:,2)./V1(:,1);
X(:,3) = V1(:,3);

N1 = sqrt(sum(V1(:,1:2).^2,2));
Nx = sqrt(sum(X(:,1:2).^2,2));
norm_pct = N1./Nx;

% Update norm = distance ratio OM/OX
V1(:,1:2) = V1(:,1:2).*norm_pct./N1;

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

% Add lower and upper disks
[Vlc,Tlc] = mesh_disk(r,edg_nb_smpl);
Vlc(:,3) = Vlc(:,3) - 0.5*h;

Vuc = Vlc;
Vuc(:,3) = Vuc(:,3) + h;
Tuc = Tlc + size(Vlc,1);

T = cat(1,T,Tlc+size(V,1),Tuc+size(V,1));
V = cat(1,V,Vlc,Vuc);

% Duplicata removal
[V,~,n] = uniquetol(V,'ByRows',true);
T = n(T);


end % mesh_cylinder