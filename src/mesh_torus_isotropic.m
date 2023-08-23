function [V, T] = mesh_torus_isotropic(R, r, nb_samples)
%
% Author & support : nicolas.douillet (at) free.fr, 2023.
%
%
% - R > r > 0 : large torus radius.
% - r > 0 : small torus radius.


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
        T(row_idx,:)   = [i*nb_samples+j i*nb_samples+j+1 (i-1)*nb_samples+j+1];
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
[V,~,n] = uniquetol(V,1e3*eps,'ByRows',true);
T = n(T);

% Remove duplicated triangles
T_sort = sort(T,2);
[~,idx,~] = unique(T_sort,'rows','stable');
T = T(idx,:);


end % mesh_torus_isotropic