function [V, T] = mesh_sphere_isotropic(r, nb_samples)
%% mesh_sphere_isotropic : function isotropically mesh a sphere.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2023-2025.
%
%
%%% Input arguments
%
% - r > 0 :      positive real scalar double, the sphere radius. Mandatory.
%
% - nb_samples : positive integer scalar double , the number of samples. Mandatory.
%
%
%%% Output arguments
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
%%% About / other information
%
% Sphere is centered on the origin, [0 0 0].
% Triangles / normals are coherently oriented and facing outward.


%% Body
u = linspace(0,pi,0.5*nb_samples)';
v = linspace(0,2*pi,nb_samples);

% Surface sampling
X = r*sin(u)*cos(v);
Y = r*sin(u)*sin(v);
Z = r*repmat(cos(u),[1,nb_samples]);

% Triplet indices for mesh facets
T = zeros(0.5*nb_samples^2,3);
row_id = 1;
i = 1;

while i < nb_samples
    
    j = 1;
    
    while j < 0.5*nb_samples      
        
        T(row_id,:)   = [(i-1)*0.5*nb_samples+j (i-1)*0.5*nb_samples+j+1 i*0.5*nb_samples+j];
        row_id = row_id + 1;
        T(row_id,:)   = [(i-1)*0.5*nb_samples+j+1 i*0.5*nb_samples+j+1 i*0.5*nb_samples+j];
        row_id = row_id + 1;
        
        j = j + 1;
        
    end
    
    % begin-end loop junction
    T(row_id,:)   = [(i-1)*0.5*nb_samples+j (i-1)*0.5*nb_samples+1 i*0.5*nb_samples+j];
    row_id = row_id + 1;
    T(row_id,:)   = [i*0.5*nb_samples+j i*0.5*nb_samples (i-1)*0.5*nb_samples+1];
    row_id = row_id + 1; 
    
    i = i + 1;
    
end

X = X(:);
Y = Y(:);
Z = Z(:);

V = [X Y Z];

% Remove duplicated vertices
tol = 1e3*eps;
[V,T] = remove_duplicated_vertices(V,T,tol);

% Remove duplicated triangles
T = remove_duplicated_triangles(T);


end % mesh_sphere_isotropic