% Regression tests: imported IDs, arbitrary coordinates, cavities, input
% validation, affine patch and end-to-end file-based DVC (MAT and CSV).
function report = test_userMesh_tetrahedron
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'func'),fullfile(root,'src'),fullfile(root,'PlotFiles'));
folder=fullfile(root,'sample_userMesh_tetrahedron');
matFile=fullfile(folder,'mesh_userMesh_tetrahedron.mat');
csvSource=struct('nodesFile',fullfile(folder,'nodes_userMesh_tetrahedron.csv'), ...
    'elementsFile',fullfile(folder,'elements_userMesh_tetrahedron.csv'));
mesh=ReadMesh_userMesh_tetrahedron(matFile); meshCSV=ReadMesh_userMesh_tetrahedron(csvSource);
assert(max(abs(mesh.coordinatesFEM-meshCSV.coordinatesFEM),[],'all')<1e-10);
assert(isequal(mesh.elementsFEM,meshCSV.elementsFEM) && any(mesh.orientationCorrected));
raw=load(matFile);
assert(isequal(mesh.inputNodeIDs,raw.nodeIDs) && isequal(mesh.inputElementIDs,raw.elementIDs));
assert(isequal(sort(mesh.inputNodeIDs(mesh.elementsFEM),2),sort(raw.elementsFEM,2)));
assert(any(mod(mesh.coordinatesFEM(:),1)~=0));
[~,V]=funTetGeometry_userMesh_tetrahedron(mesh);
assert(all(V>0));
%% Hole and interpolation: respect supplied topology, not its convex hull
[q,~,~]=funMeshVoxels_userMesh_tetrahedron(mesh,[40 40 40]);
assert(~any(all(q==[20 20 20],2)) && size(unique(q,'rows'),1)==size(q,1));
TR=triangulation(mesh.elementsFEM,mesh.coordinatesFEM);
assert(isnan(pointLocation(TR,[20 20 20])));
H=[.012 .021 -.009;-.015 -.008 .017;.006 -.011 .019];
u=mesh.coordinatesFEM*H'+[.2 -.1 .3]; U=reshape(u',[],1);
[F,Fe]=funGlobal_NodalStrainAvg3_userMesh_tetrahedron(mesh,U,1);
assert(max(abs(F-repmat(H(:),size(u,1),1)))<1e-12);
assert(max(abs(Fe-repmat(H(:)',size(Fe,1),1)),[],'all')<1e-12);
%% File-based DVC with zero initialization, no normalization for exact truth
ref=fullfile(folder,'vol_reference_tetrahedron.mat'); def=fullfile(folder,'vol_deformed_tetrahedron.mat');
p=struct('alpha',0.1,'tol',1e-4,'maxIter',60,'plotResults',false,'normalizeImages',false);
identity=main_FE_GlobalDVC_userMesh_tetrahedron(ref,ref,matFile,p);
assert(identity.info.converged && norm(identity.U)<1e-10);
r=main_FE_GlobalDVC_userMesh_tetrahedron(ref,def,matFile,p);
rc=main_FE_GlobalDVC_userMesh_tetrahedron(ref,def,csvSource,p);
assert(r.info.converged && rc.info.converged && max(abs(r.U-rc.U))<1e-8);
rmse=sqrt(mean((reshape(r.U,3,[])-[.35;-.25;.2]).^2,'all'));
assert(rmse<0.02 && all(diff(r.info.objective)<=1e-8));
assert(isequal(r.nodalDisplacement(:,1),raw.nodeIDs));
%% Check malformed imports using a temporary file owned by this test
tempFile=[tempname(root) '_userMesh_tetrahedron.mat']; cleanup=onCleanup(@() delete(tempFile));
bad=raw; bad.elementsFEM(1,1)=max(raw.nodeIDs)+1000;
save(tempFile,'-struct','bad'); expectError(@() ReadMesh_userMesh_tetrahedron(tempFile),'tetrahedron:MissingNode');
bad=raw; bad.elementsFEM(1,2)=bad.elementsFEM(1,1);
save(tempFile,'-struct','bad'); expectError(@() ReadMesh_userMesh_tetrahedron(tempFile),'tetrahedron:RepeatedNode');
bad=raw; bad.elementsFEM(end+1,:)=bad.elementsFEM(1,:); bad.elementIDs(end+1)=max(bad.elementIDs)+1;
save(tempFile,'-struct','bad'); expectError(@() ReadMesh_userMesh_tetrahedron(tempFile),'tetrahedron:DuplicateElements');
bad=raw; bad.coordinatesFEM=bad.coordinatesFEM*2+10;
save(tempFile,'-struct','bad'); mapped=ReadMesh_userMesh_tetrahedron(tempFile,[.5 .5 .5],[-5 -5 -5]);
assert(max(abs(mapped.coordinatesFEM-mesh.coordinatesFEM),[],'all')<1e-12);
outside=mesh; outside.coordinatesFEM(:,1)=outside.coordinatesFEM(:,1)+40;
expectError(@() funMeshVoxels_userMesh_tetrahedron(outside,[40 40 40]),'tetrahedron:ROI');
report=struct('passed',true,'nodes',size(u,1),'tetrahedra',numel(V), ...
    'uniqueVoxels',size(q,1),'translationRMSE',rmse,'cavityPreserved',true);
disp(report);
end
function expectError(f,id)
try
    f();
catch err
    assert(strcmp(err.identifier,id),'Unexpected validation error: %s',err.identifier); return;
end
error('Expected validation error was not raised: %s',id);
end
