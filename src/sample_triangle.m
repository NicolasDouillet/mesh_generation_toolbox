function [V, T] = sample_triangle(V1, V2, V3, nbstep, option_random, nb_points)
% sample_triangle : function to generate a triangular sampling grid of a given triangle
% following its edge directions. An element of the resulting mesh is an
% homothetic smaller version of (V1, V2, V3) triangle. Works in any
% dimension Ndim >= 2.
%
% Author & support nicolas.douillet (at) free.fr, 2017-2021.
%
%% Syntax
% [V T] = sample_triangle(V1, V2, V3)
% [V T] = sample_triangle(V1, V2, V3, nbstep)
% V = sample_triangle(V1, V2, V3, nbstep, option_random)
% V = sample_triangle(V1, V2, V3, nbstep, option_random, nb_points)
% 
%% Description
% [V T] = sample_triangle(V1, V2, V3) generates and returns a set of points
% which samples the triangle (V1, V2, V3), such that there is 20 samples
% steps in each vector direction, (V1V2) and (V1V3).
%
% [V T] = sample_triangle(V1, V2, V3, nbstep) uses nbstep steps in each
% direction, (V1V2) and (V1V3), therefore size(V,1) is always a triangular
% number.
%
% V = sample_triangle(V1, V2, V3, nbstep, true) randomly generates
% the coordinates of 200 sampling points and the resulting grid is not
% regular anymore. Valuable precision : nbstep is not taken into account in
% this case. V1, V2, and V3 are part of the set.
% 
% V = sample_triangle(V1, V2, V3, nbstep, true, nb_points) allows
% to tune the wished number of points, nb_points. V1, V2, and V3 are part
% of the set.
%
%% See also : MESHGRID, TRIMESH, LINSPACE, MESH
%
%% Input arguments
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
% - nbstep : numeric scalar, integer, nbstep > 1. The number of sampling
%            steps.
%
% - option_random : either logical *false/true or numeric, *0/1. If false / 0, output T is irrelevant.
%
% - nb_points : numeric scalar, integer. The number of points wished in
%               case option_random = true / 1.
%
%% Output arguments
%
%        | | |
% - V = [X Y Z] in dimension 3, numeric matrix, the data matrix of the sampling points coordinates. size(V) = [nb_points,Ndim].
%        | | |        
%
%         |  |  |
% - T = [i0 i1 i2] in dimension 3. numeric matrix, triangles index matrix. size(T) = [nbstep^2,3]. Is not relevant in case option_random = true / 1.
%         |  |  |
%
%% 3D Example
%
% 3D regular sampling + mesh
% V1 = [-2 3 7]';
% V2 = [8 -1 5]';
% V3 = [3 1 -3]';
% nstep = 16;
% [V T] = sample_triangle(V1, V2, V3, nstep);
% TRI = triangulation(T, V(:,1), V(:,2), V(:,3));
% figure;
% plot3(V(:,1), V(:,2), V(:,3), 'b.', 'Linewidth', 2), hold on;
% trimesh(TRI);
% colormap([0 0 1]);
% axis equal, axis tight;
%
%% 4D example
%
% Random sampling and its 3D projection
% V1 = [-2 3 7 1]';
% V2 = [8 -1 5 2]';
% V3 = [3 1 -3 -2]';
% nstep = 16;
% option_random = true;
% V = sample_triangle(V1, V2, V3, nstep, option_random, 600);
% figure;
% plot3(V(:,1), V(:,2), V(:,3), 'b.', 'Linewidth', 2), hold on;
% axis equal, axis tight;


%% Inputs parsing
assert(nargin > 2, 'Error : not enough input arguments.');
assert(nargin < 7, 'Too many input arguments.');

if nargin < 6
    nb_points = 200;
    if nargin < 5
        option_random = false;
        if nargin < 4
            nbstep = 20;
        else
            assert(isreal(nbstep) && rem(nbstep,1) == 0 && nbstep > 0, 'nbstep must be a positive integer greater or equal to 1.')
        end
    end
end


%% Create sampling grid
global Ndim;

Ndim = size(V1,1);

% (V1V2, V1V3) base
u = (V2 - V1);
v = (V3 - V1);

V = zeros(sum(1:nbstep+1), Ndim);

if ~option_random
    
    nu = u / norm(u);
    nv = v / norm(v);
    stepu = norm(u) / nbstep;
    stepv = norm(v) / nbstep;
    k = 1;
    
    % Sampling & vertices generation    
    for m = 0:nbstep
        
        for n = 0:nbstep
            
            if (m+n <= nbstep) % in (V1,V2,V3) triangle conditions ; indices # nb segments
                
                % translation vector
                tv = m*stepu*nu + n*stepv*nv;                
                V(k,:) = (V1 + tv)';
                k = k+1;
                
            end
            
        end
        
    end
    
else
    
    V = zeros(nb_points,Ndim);
    
    rand_coeff_vect = rand(2,nb_points-3);
    f = sum(rand_coeff_vect, 1) > 1;
    rand_coeff_vect(:,f) = 1 - rand_coeff_vect(:,f);
    
    % Translation vectors
    TV = repmat(rand_coeff_vect(1,:),[Ndim 1]).*repmat(u,[1 size(rand_coeff_vect,2)]) + ...
         repmat(rand_coeff_vect(2,:),[Ndim 1]).*repmat(v,[1 size(rand_coeff_vect,2)]);        
    
    for k = 1:size(TV,2)
        
        % translation vector
        tv = TV(:,k);                
        V(k,:) = (V1 + tv)';
    end
    
    % Add the triangle vertices
    V(end-2,:) = V1';
    V(end-1,:) = V2';
    V(end,:)   = V3';
    V = unique(V','rows','stable')'; % stable order for triangles indexing in the following
end


%% Index triplets list construction
T = zeros(nbstep^2,3);
row_length = 1 + nbstep;
cum_row_length = row_length;
row_idx = 1;
p = 1;

while (p <= nbstep^2 && row_length > 1)
    
     i = p;
    
    if (p < 2) % "right" triangle serie only
        
        while (i < cum_row_length)
            T(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;
            i = i +1;
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    else
        
        % Since p >= 2
        while (i < cum_row_length) % both triangle series
            
            T(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;            
            T(row_idx,:) = [i i+1 i-row_length]; % + upside-down triangles serie
            row_idx = row_idx + 1;
            
            i = i +1;
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    end
    
end

T = sort(T,2);
T = unique(T,'rows','stable');


end % sample_triangle