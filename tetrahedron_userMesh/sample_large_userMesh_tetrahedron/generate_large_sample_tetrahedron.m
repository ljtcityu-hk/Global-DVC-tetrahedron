% =========================================================
% Offline reproducible sample preparation, NOT part of the DVC workflow.
% Produces 128^3 volumes and a supplied perturbed Tet4 mesh file.
% The reading demo only loads the files already provided in this folder.
% =========================================================
function generate_large_sample_tetrahedron
folder=fileparts(mfilename('fullpath'));
sampleInfo=struct('imageSize',[128 128 128],'translation',[0.8 -0.6 0.4], ...
    'H',[.003 .002 -.001;-.0015 -.002 .001;.0008 -.0012 .0025], ...
    'centre',[64 64 64],'coordinateOrder','Img(x,y,z)','units','voxels');
%% Reference and analytically warped images, quantized to uint16
[x,y,z]=ndgrid(1:128,1:128,1:128);
texture=@(a,b,c) sin(.27*a+.13*b-.09*c)+cos(.11*a-.25*b+.17*c)+ ...
    sin(-.19*a+.07*b+.29*c)+.5*cos(.37*a+.23*b+.15*c)+ ...
    .35*sin(.17*a-.31*b-.21*c)+.25*cos(.09*a+.41*b-.27*c);
vol={uint16(30000+5000*texture(x,y,z))};
save(fullfile(folder,'vol_reference_large_tetrahedron.mat'),'vol','sampleInfo','-v7');
t=sampleInfo.translation;
vol={uint16(30000+5000*texture(x-t(1),y-t(2),z-t(3)))};
save(fullfile(folder,'vol_translation_large_tetrahedron.mat'),'vol','sampleInfo','-v7');
% y = centre + (I+H)*(x-centre) + translation, g(y)=f(x).
A=eye(3)+sampleInfo.H; centre=sampleInfo.centre';
q=(A\([x(:),y(:),z(:)]'-centre-t'))'+centre';
vol={reshape(uint16(30000+5000*texture(q(:,1),q(:,2),q(:,3))),size(x))};
save(fullfile(folder,'vol_affine_large_tetrahedron.mat'),'vol','sampleInfo','-v7');
clear x y z q vol;
%% Prepare a conforming full solid mesh for the supplied sample
[x,y,z]=ndgrid(16:12:112,16:12:112,16:12:112);
X=[x(:),y(:),z(:)]; n=size(x,1); hex=zeros((n-1)^3,8); e=0;
for k=1:n-1
    for j=1:n-1
        for i=1:n-1
            e=e+1; a=i+(j-1)*n+(k-1)*n*n;
            hex(e,:)=[a a+1 a+n+1 a+n a+n*n a+n*n+1 a+n*n+n+1 a+n*n+n];
        end
    end
end
split=[1 2 3 7;1 3 4 7;1 4 8 7;1 8 5 7;1 5 6 7;1 6 2 7];
E=zeros(6*size(hex,1),4);
for k=1:6, E(k:6:end,:)=hex(:,split(k,:)); end
rng(20261005);
interior=all(X>16 & X<112,2);
X(interior,:)=X(interior,:)+1.2*(2*rand(nnz(interior),3)-1);
% Deliberately irregular coordinates, shuffled rows and non-contiguous IDs.
order=randperm(size(X,1)); coordinatesFEM=X(order,:);
oldIDs=(10001+17*(0:size(X,1)-1))'; nodeIDs=oldIDs(order); nodeIDs=nodeIDs(:);
elementsFEM=oldIDs(E); elementIDs=(50001+11*(0:size(E,1)-1))';
elementsFEM(1:2:end,[2 3])=elementsFEM(1:2:end,[3 2]);
save(fullfile(folder,'mesh_large_userMesh_tetrahedron.mat'), ...
    'coordinatesFEM','elementsFEM','nodeIDs','elementIDs','sampleInfo');
fprintf('Prepared 128^3 volumes, %d nodes, %d tetrahedra.\n',size(X,1),size(E,1));
end
