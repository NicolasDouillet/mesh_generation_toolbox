function [V, T] = meshed_ovoid(nb_samples, option_display)
%% meshed_ovoid : function to compute, display, and save a meshed ovoid.
%
% Author & support : nicolas.douillet (at) free.fr, 2021.
%
%
% Syntax
%
% meshed_ovoid;
% meshed_ovoid(nb_samples);
% meshed_ovoid(nb_samples, option_display);
% [V,T] = meshed_ovoid(nb_samples, option_display);
%
% Description
%
% meshed_ovoid compute and display a meshed ovoid made of 64 samples
% in longitude angle and 64 samples in latitude angle. The radius of
% the bottom half sphere equals 1.
%
% meshed_ovoid(nb_samples) uses nb_samples for a half pi angle.
%
% meshed_ovoid(nb_samples, option_display) displays the result
% when option_display is set to logical *true/*1 and doesn't when it is
% set to logical false/0.
%
% [V,T] = meshed_ovoid(nb_samples, option_display) stores the resulting
% vertices coordinates in the array V, and the corresponding triplet
% indices list in the array T.
%
%
% See also SPHERE, ELLIPSOID
%
%
% Input arguments
%
% - nb_samples : positive integer double, nb_samples > 2, the number of samples for 0.5*pi angle.
% 
% - option_display : either logical *true/false or numeric 1/0.
%
%
% Output arguments
%
%       [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the vertex coordinates. Size(V) = [nb_vertices,3].
%       [ |  |  |]
%
%       [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%       [ |  |  |]


%% Input parsing
assert(nargin < 3,'Too many input arguments.');

if nargin < 2
    
    option_display = true;
    
    if nargin < 1
        
        nb_samples = 32;
        
    else
        
        assert(isnumeric(nb_samples) && nb_samples == floor(nb_samples) && nb_samples > 2,'nb_samples value must be numeric positive integer greater or equal to 3.');
        
    end
    
else
    
    assert(islogical(option_display) || isnumeric(option_display),'option_display parameter type must be either logical or numeric.');
    
end


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
[V,~,n] = unique(V,'rows','stable');
T = n(T);

% Display
if option_display
    
    disp_ovoid(V,T);
    
end


end % meshed_ovoid


%% build_triangulation subfunction
function [T] = build_triangulation(S1, S2)


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


%% disp_ovoid subfunction
function [] = disp_ovoid(V, T)


f0 = [0 0 0];
f1 = [1 1 1];

h = figure(1);
% set(h,'Position',get(0,'ScreenSize'));

trimesh(T,V(:,1),V(:,2),V(:,3),'LineWidth',2,'FaceColor','none');

ax = gca;
set(ax,'Color',f0,'XColor',f1,'YColor',f1,'ZColor',f1);
ax.Clipping = 'off';
set(gcf,'Color',f0);

colormap([0 1 0]);
view(3);
axis equal, axis tight;


end % dispovoid