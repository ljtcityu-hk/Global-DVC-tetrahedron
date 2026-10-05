% =========================================================
% Build voxel interpolation P and exact Tet4 gradient regularizer R.
% pointLocation assigns each integer voxel to ONE tetrahedron, including
% shared faces. Image term uses unit voxel weights; R integrates volume.
% This routine consumes coordinates/connectivity, not a mesh generator.
% =========================================================
function cache = funPrepareICGN3_userMesh_tetrahedron(DVCmesh,Df,imgSize)
X = DVCmesh.coordinatesFEM; E = DVCmesh.elementsFEM;
[G,V] = funTetGeometry_userMesh_tetrahedron(DVCmesh);
[q,id,N] = funMeshVoxels_userMesh_tetrahedron(DVCmesh,imgSize);
nq = size(q,1); nn = size(X,1); nodes = E(id,:);
P = sparse(repmat((1:nq)',4,1),nodes(:),N(:),nq,nn);
offset = Df.DfAxis([1 3 5]); p = q-offset;
idx = sub2ind(size(Df.DfDx),p(:,1),p(:,2),p(:,3));
df = [Df.DfDx(idx),Df.DfDy(idx),Df.DfDz(idx)];
[ri,ci,pv] = find(P);
H = sparse([ri;ri;ri],[3*ci-2;3*ci-1;3*ci], ...
    [pv.*df(ri,1);pv.*df(ri,2);pv.*df(ri,3)],nq,3*nn);
ii = zeros(16*size(E,1),1); jj = ii; vv = ii;
for e = 1:size(E,1)
    [a,b] = ndgrid(E(e,:),E(e,:)); range = 16*(e-1)+(1:16);
    Ke = V(e)*(G(:,:,e)'*G(:,:,e));
    ii(range) = a(:); jj(range) = b(:); vv(range) = Ke(:);
end
R = kron(sparse(ii,jj,vv,nn,nn),speye(3));
cache = struct('P',P,'H',H,'R',R,'points',q,'owners',id,'shape',N, ...
    'volume',V,'sampleCounts',accumarray(id,1,[size(E,1),1]), ...
    'referenceIndex',sub2ind(imgSize,q(:,1),q(:,2),q(:,3)));
end
