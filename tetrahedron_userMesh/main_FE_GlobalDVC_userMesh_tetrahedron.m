% =========================================================
% User-supplied Tet4 connectivity and arbitrary node positions -> Global DVC.
% No mesh generation is performed in this entry point or its solver.
% meshSource: MAT filename OR struct('nodesFile',...,'elementsFile',...).
% U0: 3*N column, input node-row order, voxel displacement units; default 0.
% =========================================================
function result = main_FE_GlobalDVC_userMesh_tetrahedron(referenceFile,deformedFile,meshSource,DVCpara,U0)
%% Section 1: Paths, file inputs and parameters
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'func'),fullfile(root,'src'),fullfile(root,'PlotFiles'));
if nargin<4, DVCpara=struct; end
defaults=struct('alpha',0.1,'tol',1e-4,'maxIter',100,'StrainType',0, ...
    'plotResults',true,'normalizeImages',true,'meshScale',[1 1 1],'meshOffset',[0 0 0]);
names=fieldnames(defaults);
for k=1:numel(names)
    if ~isfield(DVCpara,names{k}), DVCpara.(names{k})=defaults.(names{k}); end
end
Img1=ReadVolume_userMesh_tetrahedron(referenceFile);
Img2=ReadVolume_userMesh_tetrahedron(deformedFile);
assert(isequal(size(Img1),size(Img2)),'tetrahedron:ImageSize','Image sizes must match.');
%% Section 2: Load existing mesh and prepare image data
DVCmesh=ReadMesh_userMesh_tetrahedron(meshSource,DVCpara.meshScale,DVCpara.meshOffset);
if nargin<5 || isempty(U0), U0=zeros(3*size(DVCmesh.coordinatesFEM,1),1); end
[q,~,~]=funMeshVoxels_userMesh_tetrahedron(DVCmesh,size(Img1));
if DVCpara.normalizeImages
    % Statistics use the actual mesh domain, excluding holes and exterior.
    idx=sub2ind(size(Img1),q(:,1),q(:,2),q(:,3));
    a=Img1(idx); b=Img2(idx);
    assert(std(a)>0 && std(b)>0,'tetrahedron:Texture','Constant grayscale in mesh domain.');
    Img1=(Img1-mean(a))/std(a); Img2=(Img2-mean(b))/std(b);
end
Df=funImgGradient3_userMesh_tetrahedron(Img1,'stencil7');
%% Section 3: Global DVC iterations
[U,normOfW,timeICGN,info]=funGlobalICGN3_userMesh_tetrahedron( ...
    DVCmesh,Df,Img1,Img2,U0,DVCpara.alpha,DVCpara.tol,DVCpara.maxIter);
%% Section 4: Element gradients, nodal recovery and strain
[F,elementGradient,centres]=funGlobal_NodalStrainAvg3_userMesh_tetrahedron(DVCmesh,U,1);
strain=ComputeStrain3_userMesh_tetrahedron(F,DVCpara.StrainType);
elementStrain=reshape(ComputeStrain3_userMesh_tetrahedron(reshape(elementGradient',[],1),DVCpara.StrainType),9,[])';
result=struct('DVCmesh',DVCmesh,'DVCpara',DVCpara,'U',U,'F',F,'Strain',strain, ...
    'elementGradient',elementGradient,'elementStrain',elementStrain, ...
    'elementCentres',centres,'normOfW',normOfW,'timeICGN',timeICGN,'info',info, ...
    'referenceFile',referenceFile,'deformedFile',deformedFile, ...
    'nodalDisplacement',[DVCmesh.inputNodeIDs,reshape(U,3,[])']);
%% Section 5: Plot the supplied mesh boundary
if DVCpara.plotResults, PlotResults3_userMesh_tetrahedron(DVCmesh,U,strain); end
end
