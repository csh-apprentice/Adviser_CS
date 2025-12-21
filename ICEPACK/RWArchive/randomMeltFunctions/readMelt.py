'''
 This code snippet demonstrates how melt functions are called.
It will not run without additional code to read thickness, h, and a floating mask, floating, which are both icepack variables
and a firedrake/icepack function space, Q.
Note icepack/firedrake must be installed.
'''
# Download git repository and add to python path
import modelfunc as mf
meltLevelFig1 = 100  # Total melt in Gt/yr
meltParams = mf.inputMeltParams(f'randomMelt.{meltLevelFig1}GT.yaml') # Read appropriate melt file
#
meltProfiles = []
#
degree = 1
mesh, Q, V, opts = mf.setupMesh(f'../meshes/PigFull2017GeomFull.exp', degree=degree, meshOversample=2)
paramsFile = '../inversionInputs/PigGeometry2017-BM2.yaml'
zb, s, h, floating, grounded = mf.getModelGeometry(paramsFile, Q, smooth=True,alpha=200)
# Using the first 30 sets of parameters (works with up to 100), generate 30 random
for mr in range(0, 30):
    # Generate a melt realization, melt, on the function space Q
    print(mr)
    melt = mf.piecewiseWithDepth(h, floating, meltParams[f'random_{mr}'], Q)

    
