% =========================================================
% Global Tet4 DVC: fixed reference-gradient Hessian, additive updates.
% Image data are summed over uniquely assigned integer voxels.
% Regularization is alpha * integral |grad(u)|^2 dV.
% Same leading input/output interface as funGlobalICGN3.
% =========================================================
function [U,normOfW,timeICGN,info] = funGlobalICGN3_userMesh_tetrahedron(DVCmesh,Df,Img1,Img2,U,alpha,tol,maxIter)
if nargin<8, maxIter=100; end
validateattributes(alpha,{'numeric'},{'scalar','finite','nonnegative'});
validateattributes(tol,{'numeric'},{'scalar','finite','positive'});
validateattributes(maxIter,{'numeric'},{'scalar','integer','positive'});
assert(isequal(size(Img1),size(Img2)),'tetrahedron:Image','Image sizes must match.');
validateattributes(U,{'numeric'},{'vector','finite','real','numel',3*size(DVCmesh.coordinatesFEM,1)});
U = double(U(:)); Img1 = double(Img1); Img2 = double(Img2);
cache = funPrepareICGN3_userMesh_tetrahedron(DVCmesh,Df,size(Img1));
P = cache.P; H = cache.H; R = cache.R;
f = Img1(cache.referenceIndex);
A = H'*H + alpha*R; A = (A+A')/2;
[L,flag] = chol(A,'lower');
assert(flag==0,'tetrahedron:Singular','DVC system is not positive definite; check texture, mesh spacing and alpha.');
normOfW = zeros(maxIter,1); timeICGN = normOfW;
info = struct('converged',false,'status','maxIter','objective',zeros(maxIter+1,1), ...
    'sampleCount',size(P,1),'elementCount',size(DVCmesh.elementsFEM,1), ...
    'elementSampleCounts',cache.sampleCounts);
if any(cache.sampleCounts==0)
    warning('tetrahedron:UnsampledElements','%d elements contain no assigned integer samples; their solution relies on neighboring elements and regularization.',nnz(cache.sampleCounts==0));
end
[r,energy] = residual(U);
assert(isfinite(energy),'tetrahedron:InitialGuess','Initial displacement maps samples outside the deformed image.');
info.objective(1)=energy;
for it = 1:maxIter
    tic;
    b = H'*r-alpha*R*U;
    W = L'\(L\b);
    assert(all(isfinite(W)),'tetrahedron:Update','Nonfinite displacement update.');
    scale = 1; accepted = false;
    for ls = 1:16
        [rNew,energyNew] = residual(U+scale*W);
        if isfinite(energyNew) && energyNew <= energy+1e-12*max(1,energy)
            accepted = true; break;
        end
        scale = scale/2;
    end
    timeICGN(it)=toc;
    if ~accepted
        normOfW(it)=norm(W)/sqrt(numel(W));
        info.status='lineSearchFailed'; info.objective(it+1)=energy; break;
    end
    U=U+scale*W; r=rNew; energy=energyNew;
    normOfW(it)=norm(scale*W)/sqrt(numel(W)); info.objective(it+1)=energy;
    fprintf('tetrahedron: iteration %d, update %.6g, objective %.6g\n',it,normOfW(it),energy);
    % Test the full proposed step: a heavily damped step is not convergence.
    if norm(W)/sqrt(numel(W))<tol
        info.converged=true; info.status='converged'; break;
    end
end
normOfW=normOfW(1:it); timeICGN=timeICGN(1:it); info.objective=info.objective(1:it+1);
if ~info.converged, warning('tetrahedron:NotConverged','Solver stopped: %s.',info.status); end

    function [r0,e0] = residual(Utrial)
        u = P*reshape(Utrial,3,[])';
        g = funSampleImage_userMesh_tetrahedron(Img2,cache.points+u);
        r0 = f-g;
        e0 = 0.5*(r0'*r0+alpha*(Utrial'*(R*Utrial)));
    end
end
