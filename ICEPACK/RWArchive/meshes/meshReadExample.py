'''
 This script inputs a mesh create with Argus and creates and icepack mesh and scalar and vector function spaces. See https://icepack.github.io/.

The mesh is setup with a function available at: https://github.com/fastice/modelfunc
Using savegmsh=True should write an intermediate copy of the mesh in gmesh format.
Note icepack/firedrake must be installed.
'''
import modelfunc as mf
meshFile = 'PigFull2017GeomFull.exp'
degree = 1  # or 2
#
# Note the original file is at half resolution, and the full resolution used in the simulations
# is achieved with meshOversample=2
#
mesh, Q, V, opts = mf.setupMesh(meshFile, degree=degree, meshOversample=2, savegmsh=False)
