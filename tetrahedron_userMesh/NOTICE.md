# Source attribution

This directory contains a Tet4 extension of the Finite-Element-Based Global
Digital Volume Correlation (FE_Global_DVC) code distributed by Jin Yang.
The upstream copyright notice is:

> Copyright (c) 2020, Jin Yang.

Copyright notice for the modifications and additions:

> Copyright (c) 2026, Jiatai.

The distribution uses the BSD 2-Clause License in [license.txt](license.txt).
The original copyright notice, redistribution conditions, and disclaimer are
retained in that file.

The image-gradient routine `func/funImgGradient3_userMesh_tetrahedron.m` is
adapted from the upstream `funImgGradient3.m`. The global DVC workflow builds
on the upstream finite-element formulation. This extension adds linear
tetrahedral element support, user-supplied mesh-file input, unique voxel
assignment, tetrahedral gradient regularization and recovery, validation,
and synthetic examples.

These extensions are modifications to the upstream project; their presence
does not imply endorsement by the upstream author.

MATLAB and its toolboxes are external dependencies and are not included in
this distribution. The optional `ba_interp3` implementation is not bundled
in this directory; the supplied examples run with MATLAB `interpn` when that
implementation is unavailable. If separately redistributing third-party
components, retain and comply with their own license terms.
