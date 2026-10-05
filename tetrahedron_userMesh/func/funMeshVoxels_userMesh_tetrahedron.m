% Sample only the supplied mesh, retaining cavities/nonconvex boundaries.
% No Delaunay reconstruction. Each integer voxel gets one owner.
function [q,id,N] = funMeshVoxels_userMesh_tetrahedron(mesh,imgSize)
X=mesh.coordinatesFEM; E=mesh.elementsFEM;
assert(all(min(X,[],1)>=4) && all(max(X,[],1)<=imgSize-3), ...
    'tetrahedron:ROI','All mesh nodes must lie in the image gradient domain [4,size-3].');
lo=ceil(min(X,[],1)); hi=floor(max(X,[],1));
assert(all(lo<=hi),'tetrahedron:Samples','No integer voxels in mesh bounds.');
TR=triangulation(E,X);
points=cell(hi(3)-lo(3)+1,1); owners=points; weights=points;
[xx,yy]=ndgrid(lo(1):hi(1),lo(2):hi(2));
for k=1:numel(points)
    p=[xx(:),yy(:),repmat(lo(3)+k-1,numel(xx),1)];
    [e,n]=pointLocation(TR,p); inside=~isnan(e);
    points{k}=p(inside,:); owners{k}=e(inside); weights{k}=n(inside,:);
end
q=vertcat(points{:}); id=vertcat(owners{:}); N=vertcat(weights{:});
assert(~isempty(id),'tetrahedron:Samples','No integer voxels inside the supplied mesh.');
end
