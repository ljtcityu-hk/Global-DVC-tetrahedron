% Read GenerateVolMatfile-compatible grayscale volume; no axis permutation.
function Img = ReadVolume_userMesh_tetrahedron(filename)
data=load(filename,'vol');
assert(isfield(data,'vol'),'tetrahedron:ImageFormat','MAT image must contain vol.');
if iscell(data.vol)
    assert(isscalar(data.vol),'tetrahedron:ImageFormat','Each file must contain one image in vol{1}.');
    Img=data.vol{1};
else
    Img=data.vol;
end
validateattributes(Img,{'numeric'},{'real','finite','nonempty'});
assert(ndims(Img)==3 && all(size(Img)>=7),'tetrahedron:ImageFormat','Expected a 3-D grayscale image, at least 7 voxels per axis.');
Img=double(Img);
end
