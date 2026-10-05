% =========================================================
% Convert column-major displacement gradients to symmetric strain tensors.
% StrainType: 0 infinitesimal; 1 Euler-Almansi; 2 Green-Lagrange.
% Input and output: nine consecutive entries per point (same as H(:)).
% =========================================================
function strain = ComputeStrain3_userMesh_tetrahedron(F,StrainType)
if nargin<2, StrainType=0; end
assert(mod(numel(F),9)==0 && all(isfinite(F(:))),'tetrahedron:Gradient','Invalid displacement gradient.');
validateattributes(StrainType,{'numeric'},{'scalar','integer','>=',0,'<=',2});
strain=zeros(numel(F),1);
for k=1:9:numel(F)
    H=reshape(F(k:k+8),3,3); A=eye(3)+H;
    switch StrainType
        case 0, S=0.5*(H+H');
        case 1
            assert(rcond(A)>1e-12,'tetrahedron:Deformation','Singular deformation gradient.');
            S=0.5*(eye(3)-(A*A')\eye(3));
        case 2, S=0.5*(A'*A-eye(3));
    end
    strain(k:k+8)=S(:);
end
end
