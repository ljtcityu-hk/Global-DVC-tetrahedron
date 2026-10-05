% =========================================================
% Read an existing linear Tet4 mesh; never generate/retriangulate a mesh.
% MAT: coordinatesFEM (N x 3), elementsFEM (E x 4), optional nodeIDs.
%      The same fields may instead be inside a DVCmesh struct.
% Text pair: source=struct('nodesFile',...,'elementsFile',...).
%      Nodes: [x y z] or [nodeID x y z].
%      Elements: [n1 n2 n3 n4] or [elementID n1 n2 n3 n4].
% Connections refer to nodeIDs when present, otherwise 1-based row numbers.
% Coordinates in solver: X = inputX .* scale + offset (voxel coordinates).
% =========================================================
function mesh = ReadMesh_userMesh_tetrahedron(source,scale,offset)
if nargin<2, scale=[1 1 1]; end
if nargin<3, offset=[0 0 0]; end
validateattributes(scale,{'numeric'},{'vector','numel',3,'finite','real','nonzero'});
validateattributes(offset,{'numeric'},{'vector','numel',3,'finite','real'});
if ischar(source) || (isstring(source) && isscalar(source))
    [~,~,ext]=fileparts(source);
    assert(strcmpi(ext,'.mat'),'tetrahedron:MeshFormat','Use a MAT mesh, or a struct with nodesFile/elementsFile for CSV/TXT files.');
    data=load(source);
    if isfield(data,'DVCmesh'), data=data.DVCmesh; end
    assert(isstruct(data) && isscalar(data) && isfield(data,'coordinatesFEM') && isfield(data,'elementsFEM'), ...
        'tetrahedron:MeshFormat','MAT mesh requires coordinatesFEM and elementsFEM.');
    X=data.coordinatesFEM; E=data.elementsFEM;
    if isfield(data,'nodeIDs'), nodeIDs=data.nodeIDs(:); else, nodeIDs=(1:size(X,1))'; end
    if isfield(data,'elementIDs'), elementIDs=data.elementIDs(:); else, elementIDs=(1:size(E,1))'; end
elseif isstruct(source) && isfield(source,'nodesFile') && isfield(source,'elementsFile')
    X=readmatrix(source.nodesFile); E=readmatrix(source.elementsFile);
    assert(ismember(size(X,2),[3 4]) && ismember(size(E,2),[4 5]), ...
        'tetrahedron:MeshFormat','Text nodes need 3/4 columns; Tet4 elements need 4/5 columns. No headers.');
    if size(X,2)==4, nodeIDs=X(:,1); X=X(:,2:4); else, nodeIDs=(1:size(X,1))'; end
    if size(E,2)==5, elementIDs=E(:,1); E=E(:,2:5); else, elementIDs=(1:size(E,1))'; end
else
    error('tetrahedron:MeshFormat','Specify a MAT filename or nodesFile/elementsFile struct.');
end
validateattributes(X,{'numeric'},{'2d','ncols',3,'finite','real','nonempty'});
validateattributes(E,{'numeric'},{'2d','ncols',4,'integer','finite','real','nonempty'});
validateattributes(nodeIDs,{'numeric'},{'vector','numel',size(X,1),'integer','finite','real'});
validateattributes(elementIDs,{'numeric'},{'vector','numel',size(E,1),'integer','finite','real'});
assert(numel(unique(nodeIDs))==numel(nodeIDs) && numel(unique(elementIDs))==numel(elementIDs), ...
    'tetrahedron:DuplicateIDs','Node/element IDs must be unique.');
[found,rows]=ismember(E,nodeIDs);
assert(all(found(:)),'tetrahedron:MissingNode','Connectivity refers to an unknown node ID.');
originalX=double(X); X=originalX.*reshape(scale,1,3)+reshape(offset,1,3); E=double(rows);
assert(size(unique(X,'rows'),1)==size(X,1),'tetrahedron:DuplicateNodes','Coincident nodes are not supported.');
assert(all(all(diff(sort(E,2),1,2)>0)),'tetrahedron:RepeatedNode','An element repeats a node.');
assert(size(unique(sort(E,2),'rows'),1)==size(E,1),'tetrahedron:DuplicateElements','Duplicate tetrahedra.');
fixed=false(size(E,1),1);
for e=1:size(E,1)
    T=(X(E(e,2:4),:)-X(E(e,1),:))';
    assert(rcond(T)>1e-12 && det(T)~=0,'tetrahedron:Degenerate','Degenerate element ID %g.',elementIDs(e));
    if det(T)<0, E(e,[2 3])=E(e,[3 2]); fixed(e)=true; end
end
% A conforming interior face has two incident tetrahedra on opposite sides.
faces=[E(:,[1 2 3]);E(:,[1 2 4]);E(:,[1 3 4]);E(:,[2 3 4])];
opposite=[E(:,4);E(:,3);E(:,2);E(:,1)];
[~,~,faceID]=unique(sort(faces,2),'rows'); counts=accumarray(faceID,1);
assert(all(counts<=2),'tetrahedron:NonManifold','More than two tetrahedra share a face.');
[~,order]=sort(faceID); paired=find(diff(faceID(order))==0);
for k=paired'
    a=order(k); b=order(k+1); p=X(faces(a,:),:);
    normal=cross(p(2,:)-p(1,:),p(3,:)-p(1,:));
    sa=dot(normal,X(opposite(a),:)-p(1,:)); sb=dot(normal,X(opposite(b),:)-p(1,:));
    assert(sign(sa)~=sign(sb),'tetrahedron:Overlap','Tetrahedra sharing a face lie on the same side.');
end
mesh=struct('coordinatesFEM',X,'elementsFEM',E,'elementType','Tet4', ...
    'inputCoordinates',originalX,'inputNodeIDs',double(nodeIDs(:)), ...
    'inputElementIDs',double(elementIDs(:)),'orientationCorrected',fixed, ...
    'coordinateScale',reshape(scale,1,3),'coordinateOffset',reshape(offset,1,3), ...
    'source',source,'dirichlet',[],'neumann',[]);
funTetGeometry_userMesh_tetrahedron(mesh);
fprintf('Loaded mesh: %d nodes, %d Tet4, %d orientations corrected.\n',size(X,1),size(E,1),nnz(fixed));
end
