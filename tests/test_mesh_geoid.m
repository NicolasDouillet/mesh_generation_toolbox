% test mesh_geoid

clc;

addpath(genpath('../src'));
addpath('../data');


id = 3;
nb_it = 3;
sampling_mode = 'edge'; % ; % 'face'
[V,T] = mesh_geoid(id,nb_it,sampling_mode);

h = figure;
set(h,'Position',get(0,'ScreenSize'));
set(gcf,'Color',[0 0 0]);

t = trisurf(T,V(:,1),V(:,2),V(:,3)); shading faceted, hold on;
colormap([0 1 1]);        

t.EdgeColor = [1 0.5 0];
set(t,'LineWidth',2);

xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal;
ax = gca;
ax.Clipping = 'off';
alpha(0.5);

set(gca,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1]);