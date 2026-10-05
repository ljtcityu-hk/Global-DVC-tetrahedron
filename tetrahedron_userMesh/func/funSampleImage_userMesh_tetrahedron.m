% Sample image at coordinates in the repository's Img(x,y,z) convention.
% Reject out-of-image queries instead of using MEX border replication.
function values = funSampleImage_userMesh_tetrahedron(Img,q)
persistent useMex
if isempty(useMex), useMex = exist('ba_interp3','file')==3; end
s = size(Img);
if any(any(q<1 | q>s)) || any(~isfinite(q(:)))
    values = nan(size(q,1),1); return;
end
if useMex
    try
        values = ba_interp3(Img,q(:,2),q(:,1),q(:,3),'cubic');
        values = values(:); return;
    catch err
        useMex = false;
        warning('tetrahedron:MexFallback','ba_interp3 unavailable (%s); using interpn cubic.',err.identifier);
    end
end
values = interpn(Img,q(:,1),q(:,2),q(:,3),'cubic',NaN);
values = values(:);
end
