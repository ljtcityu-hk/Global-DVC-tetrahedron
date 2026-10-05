% Offline cylinder phantom and conforming Tet4 mesh preparation.
% No PDE Toolbox needed: triangulated disk -> prisms -> conforming Tet4.
% The DVC demo reads the resulting files; it does not generate a mesh.
function generate_cylinder_sample_tetrahedron
folder=fileparts(mfilename('fullpath'));
sampleInfo=struct('imageSize',[128 128 128],'centre',[64 64 64], ...
    'radius',40,'zRange',[16 112],'phantomRadius',44,'phantomZRange',[12 116], ...
    'translation',[.6 -.4 .3],'H',[-.003 .002 0;0 -.003 0;0 0 .01], ...
    'coordinateOrder','Img(x,y,z)','units','voxels');
%% Section 1: Smooth texture inside a cylindrical object, constant exterior
[x,y,z]=ndgrid(1:128,1:128,1:128);
texture=@(a,b,c) sin(.27*a+.13*b-.09*c)+cos(.11*a-.25*b+.17*c)+ ...
    sin(-.19*a+.07*b+.29*c)+.5*cos(.37*a+.23*b+.15*c)+ ...
    .35*sin(.17*a-.31*b-.21*c)+.25*cos(.09*a+.41*b-.27*c);
vol={phantom(x,y,z)};
save(fullfile(folder,'vol_reference_cylinder_tetrahedron.mat'),'vol','sampleInfo','-v7');
t=sampleInfo.translation;
vol={phantom(x-t(1),y-t(2),z-t(3))};
save(fullfile(folder,'vol_translation_cylinder_tetrahedron.mat'),'vol','sampleInfo','-v7');
A=eye(3)+sampleInfo.H; centre=sampleInfo.centre';
q=(A\([x(:),y(:),z(:)]'-centre-t'))'+centre';
vol={reshape(phantom(q(:,1),q(:,2),q(:,3)),size(x))};
save(fullfile(folder,'vol_stretch_cylinder_tetrahedron.mat'),'vol','sampleInfo','-v7');
clear x y z q vol;
%% Section 2: Polygonal approximation of the circular disk
rng(20261005); xy=[0 0];
for ring=1:4
    theta=(0:8*ring-1)'*(2*pi/(8*ring));
    points=10*ring*[cos(theta),sin(theta)];
    if ring<4, points=points+.6*(2*rand(size(points))-1); end
    xy=[xy;points]; %#ok<AGROW>
end
disk=delaunayTriangulation(xy); tri=sort(disk.ConnectivityList,2);
levels=16:12:112; n=size(xy,1);
X=[repmat(xy+[64 64],numel(levels),1),repelem(levels',n)];
E=zeros(3*size(tri,1)*(numel(levels)-1),4); cursor=0;
% Globally sorted disk IDs force the same diagonal on every prism side.
for layer=1:numel(levels)-1
    low=tri+(layer-1)*n; high=low+n; count=size(tri,1);
    E(cursor+(1:count),:)=[low high(:,3)]; cursor=cursor+count;
    E(cursor+(1:count),:)=[low(:,1:2) high(:,2:3)]; cursor=cursor+count;
    E(cursor+(1:count),:)=[low(:,1) high]; cursor=cursor+count;
end
% Perturb interior nodes only, preserving cylindrical surface and end planes.
inside=vecnorm(X(:,1:2)-[64 64],2,2)<39 & X(:,3)>16 & X(:,3)<112;
X(inside,:)=X(inside,:)+.4*(2*rand(nnz(inside),3)-1);
order=randperm(size(X,1)); coordinatesFEM=X(order,:);
oldIDs=(30001+19*(0:size(X,1)-1))'; nodeIDs=oldIDs(order); nodeIDs=nodeIDs(:);
elementsFEM=oldIDs(E); elementIDs=(90001+13*(0:size(E,1)-1))';
sampleInfo.circleSegments=32;
sampleInfo.polygonVolume=32/2*40^2*sin(2*pi/32)*96;
save(fullfile(folder,'mesh_cylinder_tetrahedron.mat'), ...
    'coordinatesFEM','elementsFEM','nodeIDs','elementIDs','sampleInfo');
writematrix([nodeIDs coordinatesFEM],fullfile(folder,'nodes_cylinder_tetrahedron.csv'));
writematrix([elementIDs elementsFEM],fullfile(folder,'elements_cylinder_tetrahedron.csv'));
fprintf('Prepared cylinder: %d nodes, %d tetrahedra.\n',size(X,1),size(E,1));

    function values=phantom(a,b,c)
        mask=(a-64).^2+(b-64).^2<=44^2 & c>=12 & c<=116;
        values=uint16(30000+5000*texture(a,b,c)); values(~mask)=uint16(2000);
    end
end
