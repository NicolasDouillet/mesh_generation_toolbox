function [V, T] = mesh_ovoid(nb_samples)
% mesh_ovoid : function to mesh and ovoid.
%
% Author & support : nicolas.douillet (at) free.fr, 2021-2023.
%
%
% Ovoid axis is (Oz).


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

    % Hügelschäffer egg equation
z = (sqrt(a^2 - d^2*sin(angle_vect).^2) + d*cos(angle_vect)).*cos(angle_vect);
x = b*sin(angle_vect);
y = zeros(1,numel(x));

U = cat(1,x,y,z);
V = U';          

for theta = angle_vect
    
    R = (Rmz(theta)*U)';
    V = cat(1,V,R);
    
end

S1 = numel(x);
S2 = 2*nb_samples+1;
T = build_triangulation(S1,S2);           
           
% Remove duplicated vertices
[V,~,n] = uniquetol(V,1e3*eps,'ByRows',true);
T = n(T);

% Remove duplicated triangles
T_sort = sort(T,2);
[~,idx,~] = unique(T_sort,'rows','stable');
T = T(idx,:);


end % mesh_ovoid


% build_triangulation subfunction
function T = build_triangulation(S1, S2)


r1 = cat(2,1,repelem(2:S1-1,2),S1);
r1 = reshape(r1,[2,S1-1])';
R1 = cat(2,r1,(1:S1-1)'+S1); % 1st triangle row indices
% size(R1,1) = S1-1

r2 = cat(2,1+S1,repelem(2+S1:2*S1-1,2),2*S1);
r2 = reshape(r2,[2,S1-1])';
R2 = cat(2,(2:S1)',fliplr(r2)); % 2nd triangle row indices

T = repmat(cat(1,R1,R2),[S2-1,1]);
% size(T) = 2*(S1-1)*(S2-1)

T = T + S1*repelem((0:S2-2)',2*(S1-1),3);


end % build_triangulation