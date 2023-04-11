function [V, T] = meshed_reuleaux_tetrahedron(sample_step, option_display)
%% meshed_reuleaux_tetrahedron : function to compute, display, and save a meshed Reuleaux tetrahedron. 
%
% Author & support : nicolas.douillet (at) free.fr, 2017-2022.
%
%
% Syntax
%
% meshed_reuleaux_tetrahedron;
% meshed_reuleaux_tetrahedron(sample_step);
% meshed_reuleaux_tetrahedron(sample_step, option_display);
% [V, T] = meshed_reuleaux_tetrahedron(sample_step, option_display);
%
%
% Description
%
% meshed_reuleaux_tetrahedron computes and displays the meshed Reuleaux
% tetrahedron included in the unit sphere, and which each
% edge is sampled in 8.
%
% meshed_reuleaux_tetrahedron(sample_step) uses sample_step steps.
%
% meshed_reuleaux_tetrahedron(sample_step, option_display)
% displays the result when option_display is logical *true/1, and doesn't when it is
% logical false/0.
%
% [V, T] = meshed_reuleaux_tetrahedron(sample_step, option_display) stores the resulting
% vertices coordinates in the array V, and the corresponding triplet indices list in the array T.
% 
%
% See also
%
% <https://fr.mathworks.com/help/matlab/ref/mesh.html?s_tid=srchtitle mesh> | 
% <https://fr.mathworks.com/help/matlab/ref/trimesh.html?searchHighlight=trimesh&s_tid=doc_srchtitle trimesh>
%
%
% Input arguments
%
% - sample_step : positive integer scalar, power of 2.
%
% - option_display : logical *true (1) / false (0).
%
%
% Output arguments
%
%     [ |  |  |]
% - V [Vx Vy Vz] : real matrix double, the point set. Size = [nb_vertices,3].
%     [ |  |  |]
%
%     [ |  |  |]
% - T [T1 T2 T3] : positive integer matrix double, the triangulation. Size = [nb_triangles,3].
%     [ |  |  |]


%% Input parsing
assert(nargin < 3,'Too many input arguments.');

if nargin > 0
    assert(isnumeric(sample_step) && sample_step == floor(sample_step) && sample_step > 0,'sample_step parameter must be a positive integer.');    
    if nargin  > 1                                
            assert(islogical(option_display) || isnumeric(option_display),'option_display parameter type must be either logical or numeric.');                    
    else                
        option_display = true;        
    end    
else
    sample_step = 8;    
    option_display = true;    
end


%% Body
%  Summits of original tetrahedron, living in the sphere S(O,1)
V1 = [0 0 1];
V2 = [2*sqrt(2)/3 0 -1/3];
V3 = [-sqrt(2)/3 sqrt(6)/3 -1/3];
V4 = [-sqrt(2)/3 -sqrt(6)/3 -1/3];

edge_length = norm(V1-V2); %  = 2*sqrt(6)/3

[V123, T] = sample_triangle(V2', V1', V3', sample_step, false, 200);
V_flat = V123; 

D123 = sqrt(sum((V123 - V4).^2,2)); % distance matrix
V123 = edge_length*(V123 - V4) ./ repmat(D123, [1 3]) + repmat(V4, [size(V123,1), 1]); % "inflated triangle" / opposite vertex V4

be = 1:sample_step+1;                       % bottom edge index vector
re = cumsum(sample_step+1:-1:0);            % right edge index vector
le = cat(2,1,1+cumsum(sample_step+1:-1:2)); % left edge index vector

M34 = [-sqrt(2)/3 0 -1/3]; % middle of [V3;V4] segment
D34 = sqrt(sum((V_flat(le,:) - M34).^2,2));
V_flat(le,:) = sqrt(2)*(V_flat(le,:) - M34) ./ repmat(D34, [1 3]) + repmat(M34, [size(V_flat(le,:),1), 1]);
V123(le,:) = V_flat(le,:);

M24 = [1/3/sqrt(2) -1/sqrt(6) -1/3]; % middle of [V2;V4] segment
D24 = sqrt(sum((V_flat(re,:) - M24).^2,2));
V_flat(re,:) = sqrt(2)*(V_flat(re,:) - M24) ./ repmat(D24, [1 3]) + repmat(M24, [size(V_flat(re,:),1), 1]);
V123(re,:) = V_flat(re,:);

M14 = [-1/3/sqrt(2) -1/sqrt(6) 1/3]; % middle of [V1;V4] segment
D14 = sqrt(sum((V_flat(be,:) - M14).^2,2));
V_flat(be,:) = sqrt(2)*(V_flat(be,:) - M14) ./ repmat(D14, [1 3]) + repmat(M14, [size(V_flat(be,:),1), 1]);
V123(be,:) = V_flat(be,:);


% Tetrahedron faces rotations
Rmy = @(theta) [cos(theta) 0 -sin(theta);
                0          1  0;
                sin(theta) 0  cos(theta)];

Rmz = @(theta) [cos(theta) -sin(theta) 0;
                sin(theta)  cos(theta) 0;
                0           0          1];
                        
V134 = (Rmz(2*pi/3)*V123')';
V142 = (Rmz(2*pi/3)*V134')';
V234 = (Rmy(-acos(-1/3))*Rmz(pi/3)*V142')';
            
V = [V123; V134; V142; V234];

% Triplet indices list 
T = [T;
     T+repmat(size(V123,1),   [size(T,1) size(T,2)]);...
     T+2*repmat(size(V123,1), [size(T,1) size(T,2)]);...
     T+3*repmat(size(V123,1), [size(T,1) size(T,2)])];
    

% Duplicated vertices removal (from the edges)
[V,T] = remove_duplicated_vertices(V,T);


TRI = triangulation(T, V(:,1), V(:,2), V(:,3));

%%  Display
if option_display
    
    figure;    
    trimesh(TRI), hold on;    
    colormap([0 0 1]);    
    axis square, axis equal, axis tight;
    
end


end % meshed_reuleaux_tetrahedron


%% sample_triangle subfunction
function [T, I] = sample_triangle(V0, V1, V2, nbstep, option_random, nb_points)


% Input parsing
assert(nargin > 2, 'Error : not enough input arguments.');
assert(nargin < 7, 'Too many input arguments.');

% Body
% Create sampling grid
global Ndim;

Ndim = size(V0,1);

% (V0V1, V0V2) base
u = (V1 - V0);
v = (V2 - V0);

T = zeros(sum(1:nbstep+1), Ndim);

if ~option_random
    
    nu = u / norm(u);
    nv = v / norm(v);
    stepu = norm(u) / nbstep;
    stepv = norm(v) / nbstep;
    k = 1;
    
    % Sampling & vertices generation    
    for m = 0:nbstep
        
        for n = 0:nbstep
            
            if (m+n <= nbstep) % in (V0,V1,V2) triangle conditions ; indices # nb segments
                
                % translation vector
                tv = m*stepu*nu + n*stepv*nv;                
                T(k,:) = (V0 + tv)';
                k = k+1;
                
            end
            
        end
        
    end
    
else
    
    T = zeros(nb_points,Ndim);
    
    rand_coeff_vect = rand(2,nb_points-3);
    f = sum(rand_coeff_vect, 1) > 1;
    rand_coeff_vect(:,f) = 1 - rand_coeff_vect(:,f);
    
    % Translation vectors
    TV = repmat(rand_coeff_vect(1,:),[Ndim 1]).*repmat(u,[1 size(rand_coeff_vect,2)]) + ...
         repmat(rand_coeff_vect(2,:),[Ndim 1]).*repmat(v,[1 size(rand_coeff_vect,2)]);        
    
    for k = 1:size(TV,2)
        
        % translation vector
        tv = TV(:,k);                
        T(k,:) = (V0 + tv)';
    end
    
    % Add the triangle vertices
    T(end-2,:) = V0';
    T(end-1,:) = V1';
    T(end,:)   = V2';
    T = unique(T','rows', 'stable')'; % stable order for triangles indexing in the following
end

% Index triplets list construction %
I = zeros(nbstep^2,3);
row_length = 1 + nbstep;
cum_row_length = row_length;
row_idx = 1;
p = 1;

while (p <= nbstep^2 && row_length > 1)
    
     i = p;
    
    if (p < 2) % "right" triangle serie only
        
        while (i < cum_row_length)
            I(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;
            i = i +1;
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    else
        
        % Since p >= 2
        while (i < cum_row_length) % both triangle series
            
            I(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;            
            I(row_idx,:) = [i i+1 i-row_length]; % + upside-down triangles serie
            row_idx = row_idx + 1;
            
            i = i +1;
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    end
    
end

I = sort(I, 2);
I = unique(I, 'rows', 'stable');


end % sample_triangle


%% remove_duplicated_vertices subfunction
function [V_out, T_out] = remove_duplicated_vertices(V_in, T_in)

tol = 1e4*eps;
[V_out,~,n] = uniquetol(V_in,tol,'ByRows',true);
T_out = n(T_in);

end % remove_duplicated_vertices