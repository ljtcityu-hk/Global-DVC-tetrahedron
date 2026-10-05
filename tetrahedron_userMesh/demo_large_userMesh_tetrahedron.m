% =========================================================
% Larger file-based benchmark: 128^3 images, 729 nodes, 3072 Tet4.
% Runs translation and affine deformation from zero initial displacement.
% This demo reads pre-existing images/mesh; it never generates a mesh.
% =========================================================
clear; clc;
root=fileparts(mfilename('fullpath')); addpath(root,fullfile(root,'PlotFiles'));
%% Section 1: Existing image and mesh files
sampleFolder=fullfile(root,'sample_large_userMesh_tetrahedron');
referenceFile=fullfile(sampleFolder,'vol_reference_large_tetrahedron.mat');
deformedFiles={fullfile(sampleFolder,'vol_translation_large_tetrahedron.mat'), ...
    fullfile(sampleFolder,'vol_affine_large_tetrahedron.mat')};
meshSource=fullfile(sampleFolder,'mesh_large_userMesh_tetrahedron.mat');
metadata=load(meshSource,'sampleInfo'); sampleInfo=metadata.sampleInfo;
caseNames={'translation','affine'};
%% Section 2: Parameters
% Same photometric scale in all supplied images, so do not normalize each
% image independently. Alpha is in this uint16 grayscale scale.
DVCpara=struct('alpha',0.1,'tol',1e-4,'maxIter',60,'StrainType',0, ...
    'normalizeImages',false,'plotResults',false, ...
    'meshScale',[1 1 1],'meshOffset',[0 0 0]);
results=cell(2,1); benchmark=struct([]);
%% Section 3: Solve both cases; truth is used ONLY after solving
for caseIndex=1:2
    fprintf('\n--- Large Tet4 case: %s ---\n',caseNames{caseIndex});
    totalTimer=tic;
    result=main_FE_GlobalDVC_userMesh_tetrahedron( ...
        referenceFile,deformedFiles{caseIndex},meshSource,DVCpara,[]);
    totalSeconds=toc(totalTimer);
    X=result.DVCmesh.coordinatesFEM;
    if caseIndex==1, H=zeros(3); else, H=sampleInfo.H; end
    truth=(X-sampleInfo.centre)*H'+sampleInfo.translation;
    recovered=reshape(result.U,3,[])'; errors=recovered-truth;
    rmse=sqrt(mean(errors.^2,'all'));
    maxVectorError=max(sqrt(sum(errors.^2,2)));
    gradientError=result.elementGradient-H(:)';
    entry=struct('caseName',caseNames{caseIndex}, ...
        'imageSize',sampleInfo.imageSize,'imageVoxelCount',prod(sampleInfo.imageSize), ...
        'sampleCount',result.info.sampleCount,'nodeCount',size(X,1), ...
        'elementCount',size(result.DVCmesh.elementsFEM,1),'dofs',numel(result.U), ...
        'converged',result.info.converged,'iterations',numel(result.normOfW), ...
        'displacementRMSE',rmse,'maxVectorError',maxVectorError, ...
        'elementGradientRMSE',sqrt(mean(gradientError.^2,'all')), ...
        'totalSeconds',totalSeconds,'iterationSeconds',sum(result.timeICGN));
    if caseIndex==1, benchmark=entry; else, benchmark(caseIndex)=entry; end
    disp(benchmark(caseIndex));
    assert(result.info.converged,'tetrahedron:Benchmark','Case did not converge.');
    assert(rmse<0.03,'tetrahedron:Benchmark','Displacement RMSE exceeds 0.03 voxels.');
    assert(result.info.sampleCount==97^3,'tetrahedron:Benchmark','Unexpected full-solid voxel coverage.');
    assert(all(result.info.elementSampleCounts>0),'tetrahedron:Benchmark','Unsampled elements.');
    results{caseIndex}=result;
    %% Section 4: Cutaway visualization and exports
    fig=PlotCutaway_userMesh_tetrahedron(result.DVCmesh,result.U,truth, ...
        ['128^3 Tet4 - ' caseNames{caseIndex}]);
    exportgraphics(fig,fullfile(root,['large_' caseNames{caseIndex} '_tetrahedron.png']),'Resolution',140);
    writematrix(result.nodalDisplacement, ...
        fullfile(root,['displacement_large_' caseNames{caseIndex} '_tetrahedron.csv']));
end
save(fullfile(root,'results_large_userMesh_tetrahedron.mat'),'results','benchmark','sampleInfo','DVCpara');
writetable(struct2table(benchmark),fullfile(root,'benchmark_large_userMesh_tetrahedron.csv'));
