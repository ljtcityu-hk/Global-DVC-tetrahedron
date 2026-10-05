% Plot the actual Tet4 boundary connectivity (no new Delaunay mesh).
function figures = PlotResults3_userMesh_tetrahedron(DVCmesh,U,strain)
TR=triangulation(DVCmesh.elementsFEM,DVCmesh.coordinatesFEM);
faces=freeBoundary(TR); X=DVCmesh.coordinatesFEM;
figures=gobjects(2,1);
figures(1)=figure('Name','tetrahedron displacement','Color','w');
labels={'u','v','w'};
for k=1:3
    subplot(1,3,k); draw(U(k:3:end)); title(labels{k});
end
figures(2)=figure('Name','tetrahedron strain','Color','w');
components=[1 5 9 4 7 8]; labels={'e11','e22','e33','e12','e13','e23'};
for k=1:6
    subplot(2,3,k); draw(strain(components(k):9:end)); title(labels{k});
end
    function draw(values)
        patch('Faces',faces,'Vertices',X,'FaceVertexCData',values, ...
            'FaceColor','interp','EdgeColor','none');
        axis equal tight; view(3); colorbar;
        xlabel('x (voxels)'); ylabel('y (voxels)'); zlabel('z (voxels)');
    end
end
