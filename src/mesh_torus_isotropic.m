function [V, T] = mesh_torus_isotropic(R, r, nb_samples)
% mesh_torus_isotropic : function isotropically mesh a torus.
%
% Author & support : nicolas.douillet (at) free.fr, 2023.
%
%
%%% Input arguments :
%
% - R > r > 0 :  positive real scalar double, the torus large radius.
% - r > 0 :      positive real scalar double, the torus small radius.
% - nb_samples : positive integer scalar, the number of samples.
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


u = linspace(0,2*pi,nb_samples)';
v = linspace(0,2*pi,nb_samples);

% Surface sampling
X = repmat(R*cos(v),[nb_samples,1]) + r*sin(u)*cos(v);
Y = repmat(R*sin(v),[nb_samples,1]) + r*sin(u)*sin(v);
Z = repmat(r*cos(u),[1,nb_samples]);

% Triplet indices for mesh facets
T = zeros(nb_samples*nb_samples, 3);
row_idx = 1;
i = 1;

while i < nb_samples 
    
    j = 1;
    
    while j < nb_samples      
        
        T(row_idx,:)   = [(i-1)*nb_samples+j (i-1)*nb_samples+j+1 i*nb_samples+j];
        row_idx = row_idx + 1;
        T(row_idx,:)   = [(i-1)*nb_samples+j+1 i*nb_samples+j+1 i*nb_samples+j];
        row_idx = row_idx + 1;
        
        j = j + 1;
        
    end
    
    % begin-end loop junction
    T(row_idx,:)   = [(i-1)*nb_samples+j (i-1)*nb_samples+1 i*nb_samples+j];
    row_idx = row_idx + 1;
    T(row_idx,:)   = [i*nb_samples+j i*nb_samples (i-1)*nb_samples+1];
    row_idx = row_idx + 1;
    
    i = i + 1;
    
end

X = X(:);
Y = Y(:);
Z = Z(:);

V = [X Y Z];

% Remove duplicated vertices
[V,T] = remove_duplicated_vertices(V,T);

% Remove duplicated triangles
T = remove_duplicated_triangles(T);


end % mesh_torus_isotropic