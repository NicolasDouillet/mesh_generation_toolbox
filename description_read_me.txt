Please first check the doc tab on the right to get some relevant and functional examples of this toobox functions.

Feel free to check and download the mesh processing toolbox ( https://fr.mathworks.com/matlabcentral/fileexchange/77004-mesh-processing-toolbox?s_tid=srchtitle ) for many further meshing relative functions.


%% DESCRIPTION

This package is a mesh generation toolbox which on the long term aim at providing a command line mesh generation lab in Matlab(R) console. It is designed to deal with and generate 3D triangular meshes.


%% HELP

A basic help is included in the header of each source file. It especially includes input and output arguments precise descriptions (role, class, size, etc).
Just like for any Matlab (R) function, typewrite "help my_mesh_generation_file" in Matlab console to access  it.


%% DATA FORMATS & HYPOTHESIS ON THE MESH

Most of the functions included take very common and widely used data structures as inputs and outputs :

- V : the vertex set / point cloud. Real matrix of double. size(V) = [nb_vertex,3].
- T : the triangulation / triangle set. Positive integer matrix of double. size(T) = [nb_triangles,3].
- E : the edge set. Positive integer matrix of double. size(T) = [nb_edges,2].

where :

- A vextex is a 1 x 3 row double vector of real numbers.
- A triangle (or a triplet) is a 1 x 3 row double vector of positive integers.
- An edge is a 1 x 2 row double vector of positive integers.

By default, and unless exceptions, vertex and triangle arguments are index based.

Another common argument is ngb_degre which corresponds to the neighborhood degre whished on the mesh, and which is used to find one vertex neighborhood.
For commun usages, this value is in the range |[1 ; 4]|. Tune it relatively to the estimated local curvature of your point set.


%% TESTS

Use .mat data files provided in /data for test and example files.
Most of the functions and every important ones have been tested in a dedicated file named : test_my_function.m.
Note that no mesh reader or writer is provided in this toolbox since there already exist enough satisfying ones coded in Matlab.
Look for : read_ply.m, write_ply.m, read_off.m, write_off.m, plyread.m, plywrite.m
Then to create your own .mat file for vertex set V and triangle set T, just use the command save('path_to_my_file/my_file.mat','V','T');


%% COPYRIGHT & SUPPORT

All the code included in this toolbox is the result of my unique personal own work and effort in 2020-2023, and going on for upcoming updates.

Each one of the algorithms / functions included have been independently tested, however I cannot provide any warrantee of any kind about them. Use them at your own risks.
Downloading and using this toolbox or just part of it assume you to have read and accept all the condition in this current description. 

This toolbox and its content is free of use and distribution with the following condition :
this description_read_me file must be included as well as each function header must be preserved.

SELLING THIS WHOLE TOOLBOX OR EVEN PARTS OF IT IS STRICLY PROHIBITED.

Modification of any kind are done under your own, only, and unique responsability.

Please report me any bug (with data set used and Matlab(R) code attached) or suggestion at nicolas.douillet (at) free.fr


%% KNOWN LIMITATIONS

Three core mesh generation functions are exploitable at the moment :

- mesh_geoid
- discrete_contour_mesh_patch
- volumic_mesh_from_convex_set_mesh

The function build_volumic_mesh_from_convex_set_mesh will of course create non manifold edges. It will also lose orientation of the face normals.


%% MISC INFORMATION

Except for the ovoid shape, all the other meshes have their normals / faces coherently oriented outward.

Most of the time, I did my best to make function names pretty explicit in english.

By default, vertex and face normals are normalized at the same time they are computed.

Basic 3D mathematical computation algorithm (like point_to_plane_distance) are also independently available with their documentations in my file exchange contributions.

I especially thank William V, Binbin Qi, for what they taught me while solving my cody problems.
Matlab users, your advices and tips to improve and speed up my algorithms are welcome !

If you can’t see the mesh while plotting it, try ‘shading flat’ (it may be shadowed by its numerous dark edges).

Since I am not native english speaker, please forgive my langage approximations.

Matlab release version used for development and tests : R2019b.


Last update : 24 / 08 / 2023.