% =========================================================
% Constant Tet4 displacement gradients and volume-weighted nodal recovery.
% F contains grad(u), NOT I+grad(u). Nine entries per node, column-major:
% [du/dx,dv/dx,dw/dx,du/dy,dv/dy,dw/dy,du/dz,dv/dz,dw/dz].
% One output integration point (centroid) per Tet4; order kept for API use.
% =========================================================
function [F,StrainGaussPt,CoordsGaussPt] = funGlobal_NodalStrainAvg3_userMesh_tetrahedron(DVCmesh,U,GaussPtOrder)
if nargin<3, GaussPtOrder=1; end
validateattributes(GaussPtOrder,{'numeric'},{'scalar','positive'});
[G,V,CoordsGaussPt] = funTetGeometry_userMesh_tetrahedron(DVCmesh);
E = DVCmesh.elementsFEM; nn = size(DVCmesh.coordinatesFEM,1);
validateattributes(U,{'numeric'},{'vector','finite','numel',3*nn});
u = reshape(U,3,[])'; sums=zeros(nn,9); mass=zeros(nn,1);
StrainGaussPt=zeros(size(E,1),9);
for e=1:size(E,1)
    H = u(E(e,:),:)'*G(:,:,e)';
    StrainGaussPt(e,:)=H(:)';
    sums(E(e,:),:) = sums(E(e,:),:)+V(e)*repmat(H(:)',4,1);
    mass(E(e,:)) = mass(E(e,:))+V(e);
end
nodal=sums./mass; F=reshape(nodal',[],1);
end
