% =========================================================
% Geometry shared by Tet4 assembly and strain recovery.
% gradN(:,:,e): 3-by-4, derivatives with respect to image x,y,z.
% Coordinates follow Img(x,y,z), in voxel units.
% =========================================================
function [gradN,volume,centres] = funTetGeometry_userMesh_tetrahedron(DVCmesh)
X = DVCmesh.coordinatesFEM; E = DVCmesh.elementsFEM;
validateattributes(X,{'double'},{'2d','ncols',3,'finite','real','nonempty'});
validateattributes(E,{'double'},{'2d','ncols',4,'integer','positive','nonempty','<=',size(X,1)});
assert(numel(unique(E(:)))==size(X,1),'tetrahedron:UnusedNodes','Every node must belong to an element.');
gradN = zeros(3,4,size(E,1)); volume = zeros(size(E,1),1); centres = zeros(size(E,1),3);
D = [-1 1 0 0;-1 0 1 0;-1 0 0 1];
for e = 1:size(E,1)
    Xe = X(E(e,:),:); T = (Xe(2:4,:)-Xe(1,:))';
    assert(det(T)>0 && rcond(T)>1e-12,'tetrahedron:Element', ...
        'Element %d has nonpositive orientation or is degenerate.',e);
    gradN(:,:,e) = T'\D;
    volume(e) = det(T)/6; centres(e,:) = mean(Xe,1);
end
end
