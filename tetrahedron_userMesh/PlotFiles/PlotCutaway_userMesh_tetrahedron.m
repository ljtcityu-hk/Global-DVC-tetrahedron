% Element cutaway: keep complete tetrahedra with centroid x <= cutX.
% This is a visualization selection only; the solve always uses the full mesh.
function fig=PlotCutaway_userMesh_tetrahedron(mesh,U,truth,titleText)
X=mesh.coordinatesFEM; E=mesh.elementsFEM;
centres=(X(E(:,1),:)+X(E(:,2),:)+X(E(:,3),:)+X(E(:,4),:))/4;
cutX=mean([min(X(:,1)),max(X(:,1))]); E=E(centres(:,1)<=cutX,:);
used=unique(E(:)); map=zeros(size(X,1),1); map(used)=1:numel(used);
faces=freeBoundary(triangulation(map(E),X(used,:)));
u=reshape(U,3,[])'; err=sqrt(sum((u-truth).^2,2));
fig=figure('Name',titleText,'Color','w','Position',[100 100 1250 520]);
subplot(1,2,1); draw(u(used,1)); title('Recovered u; element cutaway');
subplot(1,2,2); draw(err(used)); title('Displacement vector error (voxels)');
sgtitle(titleText);
    function draw(values)
        patch('Faces',faces,'Vertices',X(used,:),'FaceVertexCData',values, ...
            'FaceColor','interp','EdgeColor',[.25 .25 .25],'LineWidth',.25);
        axis equal tight; view(125,25); colorbar;
        xlabel('x (voxels)'); ylabel('y (voxels)'); zlabel('z (voxels)');
    end
end
