% =========================================================
% File-based cylinder benchmark: rigid translation and axial stretch/shear.
% Reads pre-generated cylinder mesh/volumes; both solves start from zero.
% =========================================================
clear; clc;
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'func'),fullfile(root,'src'),fullfile(root,'PlotFiles'));
%% Section 1: Read sample files
folder=fullfile(root,'sample_cylinder_tetrahedron');
referenceFile=fullfile(folder,'vol_reference_cylinder_tetrahedron.mat');
deformedFiles={fullfile(folder,'vol_translation_cylinder_tetrahedron.mat'), ...
    fullfile(folder,'vol_stretch_cylinder_tetrahedron.mat')};
meshSource=fullfile(folder,'mesh_cylinder_tetrahedron.mat');
data=load(meshSource,'sampleInfo'); sampleInfo=data.sampleInfo;
mesh=ReadMesh_userMesh_tetrahedron(meshSource);
%% Section 2: Check volume, conformity and face connectivity
X=mesh.coordinatesFEM; E=mesh.elementsFEM; TR=triangulation(E,X);
[~,volume]=funTetGeometry_userMesh_tetrahedron(mesh);
assert(abs(sum(volume)-sampleInfo.polygonVolume)/sampleInfo.polygonVolume<1e-10);
nb=neighbors(TR); [row,col]=find(~isnan(nb));
adjacent=nb(sub2ind(size(nb),row,col)); g=graph(row,adjacent,[],size(E,1));
assert(max(conncomp(g))==1,'tetrahedron:Cylinder','Disconnected volume mesh.');
faces=freeBoundary(TR);
edges=sort([faces(:,[1 2]);faces(:,[2 3]);faces(:,[3 1])],2);
[~,~,ids]=unique(edges,'rows'); assert(all(accumarray(ids,1)==2));
bg=graph(edges(:,1),edges(:,2),[],size(X,1)); bc=conncomp(bg);
assert(numel(unique(bc(unique(faces(:)))))==1,'tetrahedron:Cylinder','Unexpected boundary components.');
assert(~isnan(pointLocation(TR,[64 64 64])),'tetrahedron:Cylinder','Unexpected central cavity.');
fprintf('Cylinder connected; all element volumes positive; total volume %.3f voxels^3.\n',sum(volume));
%% Section 3: Save a surface mesh / phantom slice overview
fig=figure('Color','w','Position',[100 100 1150 500]);
subplot(1,2,1);
patch('Faces',faces,'Vertices',X,'FaceColor',[.65 .82 .95], ...
    'EdgeColor',[.2 .3 .4],'LineWidth',.35); axis equal tight; view(125,25);
xlabel('x (voxels)'); ylabel('y (voxels)'); zlabel('z (voxels)'); title('Cylinder Tet4 boundary');
img=ReadVolume_userMesh_tetrahedron(referenceFile);
subplot(1,2,2); imagesc(1:128,1:128,img(:,:,64)'); axis image xy;
xlabel('x (voxels)'); ylabel('y (voxels)'); title('Reference image, z = 64'); colormap(gca,gray); colorbar;
hold on; theta=linspace(0,2*pi,129); plot(64+40*cos(theta),64+40*sin(theta),'r-','LineWidth',1.4);
exportgraphics(fig,fullfile(root,'cylinder_mesh_tetrahedron.png'),'Resolution',140);
%% Section 4: DVC parameters and solves (no FFT or known displacement seed)
p=struct('alpha',.1,'tol',1e-4,'maxIter',60,'StrainType',2, ...
    'normalizeImages',false,'plotResults',false);
caseNames={'translation','stretch'}; results=cell(2,1); benchmark=struct([]);
for k=1:2
    timer=tic;
    result=main_FE_GlobalDVC_userMesh_tetrahedron(referenceFile,deformedFiles{k},meshSource,p,[]);
    elapsed=toc(timer);
    if k==1, H=zeros(3); else, H=sampleInfo.H; end
    truth=(X-sampleInfo.centre)*H'+sampleInfo.translation;
    errorU=reshape(result.U,3,[])'-truth;
    expectedStrain=ComputeStrain3_userMesh_tetrahedron(H(:),2);
    entry=struct('caseName',caseNames{k},'imageVoxels',prod(sampleInfo.imageSize), ...
        'sampleCount',result.info.sampleCount,'nodes',size(X,1),'tetrahedra',size(E,1), ...
        'converged',result.info.converged,'iterations',numel(result.normOfW), ...
        'displacementRMSE',sqrt(mean(errorU.^2,'all')), ...
        'maxVectorError',max(vecnorm(errorU,2,2)), ...
        'elementGradientRMSE',sqrt(mean((result.elementGradient-H(:)').^2,'all')), ...
        'elementStrainRMSE',sqrt(mean((result.elementStrain-expectedStrain').^2,'all')), ...
        'totalSeconds',elapsed);
    if k==1, benchmark=entry; else, benchmark(k)=entry; end
    disp(entry);
    assert(result.info.converged && entry.displacementRMSE<.03,'tetrahedron:Cylinder','DVC accuracy/convergence check failed.');
    assert(entry.elementGradientRMSE<.002,'tetrahedron:Cylinder','Gradient recovery check failed.');
    assert(all(result.info.elementSampleCounts>0),'tetrahedron:Cylinder','Some elements have no samples.');
    assert(all(diff(result.info.objective)<=1e-5),'tetrahedron:Cylinder','Objective increased.');
    results{k}=result;
    %% Section 5: Cutaway results and export with original node IDs
    fig=PlotCutaway_userMesh_tetrahedron(mesh,result.U,truth,['Cylinder - ' caseNames{k}]);
    exportgraphics(fig,fullfile(root,['cylinder_' caseNames{k} '_tetrahedron.png']),'Resolution',140);
    writematrix(result.nodalDisplacement,fullfile(root,['displacement_cylinder_' caseNames{k} '_tetrahedron.csv']));
end
save(fullfile(root,'results_cylinder_tetrahedron.mat'),'results','benchmark','sampleInfo');
writetable(struct2table(benchmark),fullfile(root,'benchmark_cylinder_tetrahedron.csv'));
