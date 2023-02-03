# 3D Semi-body conformal curvilinear structured grid generator for axisymmetric body
Grid points are structured as MATLAB function ndgrid() [not meshgrid()] in curvilinear coordinates (not Cartesian grid). Grid can be visualized in TecPlot.   
At first, grids are created on the Cartesian coordinate system, then a transformation is applied the conform the grid to the surface. The main idea is to adjust the radius of all the grid points on Cartesian grid while keeping the angle same.
