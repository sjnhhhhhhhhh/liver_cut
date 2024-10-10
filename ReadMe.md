#任务简介
目前完成了肝脏八段的分割，功能包含血管骨架化、以及根据分段后的肝门静脉寻找分割关键点（基于delaunay三角化和vonoroi图）进行肝脏段切割

#Introduction
This program performs liver segment segmentation using a pre-segmented hepatic portal vein. The function for locating key cutting points is based on Delaunay 3D triangulation and Voronoi diagrams, with the algorithm provided by the CGAL library. The cutting plane reconstruction uses vtkSurfaceReconstruction to fit the cutting point set.

To use this program, please ensure that you have the following libraries installed: VTK/Eigen/CGAL.

Additionally, this project provides stunning results.
