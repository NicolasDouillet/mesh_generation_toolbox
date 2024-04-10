function [V, T] = mesh_triangle(V1, V2, V3, nb_samples)
% mesh_triangle : function to sample one given (V1V2V3) triangle using nb_samples^2
% homothetic smaller versions of it. Works in any dimension Ndim >= 2.
%
% Author : nicotangente (at) free.fr, 2017-2024.
%
%
% Syntax
% [V T] = mesh_triangle(V1, V2, V3)
% [V T] = mesh_triangle(V1, V2, V3, nb_samples)
% 
% Description
% [V T] = mesh_triangle(V1, V2, V3) generates and returns a set of points
% which samples the triangle (V1, V2, V3), such that there is 20 samples
% steps in each vector direction, (V1V2) and (V1V3).
%
% [V T] = mesh_triangle(V1, V2, V3, nb_samples) uses nb_samples steps in each
% direction, (V1V2) and (V1V3), therefore size(V,1) is always a triangular
% number.
%
% See also : MESHGRID, TRIMESH, LINSPACE, MESH
%
% Input arguments
%
%        [V1x]
% - V1 = [V1y] : numeric column vector, size(V1) = [Ndim, 1]. The 1st vertex coordinates.
%        [V1z]
%
%        [V1x]
% - V2 = [V2y] : numeric column vector, size(V2) = [Ndim, 1]. The 2nd vertex coordinates.
%        [V2z]
%
%        [V3x]
% - V3 = [V3y] : numeric column vector, size(V3) = [Ndim, 1]. The 3rd vertex coordinates.
%        [V3z]
%
% - nb_samples : numeric scalar, integer, nb_samples > 1. The number of samples.
%
% Output arguments
%
%        | | |
% - V = [X Y Z] in dimension 3, numeric matrix, the data matrix of the sampling points coordinates. size(V) = [sum(1:1+nb_samples),Ndim].
%        | | |        
%
%         |  |  |
% - T = [i0 i1 i2] in dimension 3. numeric matrix, triangles index matrix. size(T) = [nb_samples^2,3].
%         |  |  |
%
% 3D Example
%
% 3D regular sampling + mesh
% V1 = [-2 3 7]';
% V2 = [8 -1 5]';
% V3 = [3 1 -3]';
% nstep = 16;
% [V T] = mesh_triangle(V1, V2, V3, nstep);
% TRI = triangulation(T, V(:,1), V(:,2), V(:,3));
% figure;
% plot3(V(:,1), V(:,2), V(:,3), 'b.', 'Linewidth', 2), hold on;
% trimesh(TRI);
% colormap([0 0 1]);
% axis equal, axis tight;
%
%
% About / other informations
%
% Triangles / normals are coherently oriented.


% Inputs parsing
assert(nargin > 2, 'Error : not enough input arguments.');
assert(nargin < 5, 'Too many input arguments.');

if nargin < 4
    nb_samples = 20;
else
    assert(isreal(nb_samples) && rem(nb_samples,1) == 0 && nb_samples > 0, 'nb_samples must be a positive integer greater or equal to 1.')
end


% Create sampling grid
global Ndim;

Ndim = size(V1,1);

% (V1V2, V1V3) base
u = (V2 - V1);
v = (V3 - V1);

V = zeros(sum(1:nb_samples+1), Ndim);

    
nu = u / norm(u);
nv = v / norm(v);
stepu = norm(u) / nb_samples;
stepv = norm(v) / nb_samples;
k = 1;

% Sampling & vertices generation
for m = 0:nb_samples
    
    for n = 0:nb_samples
        
        if m+n <= nb_samples % in (V1,V2,V3) triangle conditions ; indices # nb segments
            
            % translation vector
            tv = m*stepu*nu + n*stepv*nv;
            V(k,:) = (V1 + tv)';
            k = k+1;
            
        end
        
    end
    
end


% Index triplets list construction
T = zeros(nb_samples^2,3);
row_length = 1 + nb_samples;
cum_row_length = row_length;
row_idx = 1;
p = 1;

while p <= nb_samples^2 && row_length > 1
    
     i = p;
    
    if p < 2 % "right" triangle serie only
        
        while i < cum_row_length
            
            T(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;
            i = i +1;
            
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    else
        
        % Since p >= 2
        while i < cum_row_length % both triangle series
            
            T(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;            
            T(row_idx,:) = [i i-row_length i+1]; % + upside-down triangles serie
            row_idx = row_idx + 1;
            
            i = i +1;
            
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    end
    
end


end % mesh_triangle