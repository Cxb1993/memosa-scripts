#!/usr/bin/env python
import pdb
import sys
sys.setdlopenflags(0x100|0x2)

import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import numpy as np
import pdb

atype = 'double'
#atype = 'tangent'

if atype == 'double':
    import fvm.models_atyped_double as models
    import fvm.exporters_atyped_double as exporters
elif atype == 'tangent':
    import fvm.models_atyped_tangent_double as models
    import fvm.exporters_atyped_tangent_double as exporters


from FluentCase import FluentCase

#fvmbaseExt.enableDebug("cdtor")

fileBase = None
numIterations = 10
#fileBase = "/home/sm/app-memosa/src/fvm/test/cav32"
#fileBase = "/home/sm/a/data/wj"

def usage():
    print "Usage: %s filebase [outfilename]" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-prism.dat if it is not specified."
    sys.exit(1)

def advance(fmodel,niter):
    for i in range(0,niter):
        try:
            fmodel.advance(1)
        except KeyboardInterrupt:
            break

# change as needed

outfile = None
if __name__ == '__main__' and fileBase is None:
    if len(sys.argv) < 2:
        usage()
    fileBase = sys.argv[1]
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

if outfile == None:
    outfile = fileBase+"-prism.dat"
    
reader = FluentCase(fileBase+".cas")

#import debug
reader.read();

meshes = reader.getMeshList()

import time
t0 = time.time()

geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

if atype == 'tangent':
    metricsCalculator.setTangentCoords(0,7,1)

flowFields =  models.FlowFields('flow')

fmodel = models.FlowModelA(geomFields,flowFields,meshes)

reader.importFlowBCs(fmodel,meshes)
bcMap = fmodel.getBCMap()

bc = bcMap[7]
bc.setVar("specifiedXVelocity",0.0)

for i in [5,6,8,9,10]:
    bc = bcMap[i]
    #bc.bcType = 'Symmetry'
    bc.bcType = 'NoSlipWall'
fmodel.printBCs()

momSolver = fvmbaseExt.AMG()
momSolver.relativeTolerance = 1e-1
momSolver.nMaxIterations = 20
momSolver.maxCoarseLevels=20
momSolver.verbosity=0

contSolver = fvmbaseExt.AMG()
#pc = fvmbaseExt.AMG()
#pc.verbosity=0
#contSolver = fvmbaseExt.BCGStab()
#contSolver.preconditioner = pc
contSolver.relativeTolerance = 1e-1
contSolver.nMaxIterations = 20
contSolver.verbosity=0
contSolver.maxCoarseLevels=20

foptions = fmodel.getOptions()

foptions.momentumLinearSolver = momSolver
foptions.pressureLinearSolver = contSolver

foptions.momentumTolerance=1e-3
foptions.continuityTolerance=1e-3
foptions.setVar("momentumURF",0.7)
foptions.setVar("pressureURF",0.3)
foptions.printNormalizedResiduals=False
foptions.transient = True
foptions.setVar("timeStep",0.01)

"""
if atype=='tangent':
    vcMap = fmodel.getVCMap()
    for i,vc in vcMap.iteritems():
        print vc.getVar('viscosity')
        vc.setVar('viscosity',(1.7894e-5,1))
"""
numTimeSteps = 0
fmodel.init()

cells = meshes[0].getCells()
faces = meshes[0].getFaces()
faceCells = meshes[0].getAllFaceCells()
source = flowFields.source[cells].asNumPyArray()
rho = flowFields.density[cells].asNumPyArray()
CellCoord = geomFields.coordinate[cells].asNumPyArray()
cellPressure = flowFields.pressure[cells].asNumPyArray()
facePressure = flowFields.pressure[faces].asNumPyArray()
grav = np.array([0.0,-9.81,0.0])

source[:] = 0.0

for c in range(cells.getSelfCount()):
    source[c] = rho[c]*grav
    cellPressure[c] = rho[c]*grav[1]*CellCoord[c,1]

for f in range(faces.getSelfCount()):
    c0 = faceCells(f,0)                                                                                                  
    c1 = faceCells(f,1)
    facePressure[f] = 0.5*(cellPressure[c0]+cellPressure[c1])

pdb.set_trace()
writer = exporters.VTKWriterA(geomFields,meshes,"flow-" + str(numTimeSteps) + ".vtk",
                              "FlowField",False,0)
writer.init()
writer.writeScalarField(flowFields.pressure,"Pressure")
writer.writeVectorField(flowFields.velocity,"Velocity")
writer.finish()

while (numTimeSteps < 100):
    #advance(fmodel,numIterations)
    fmodel.advance(50)
    fmodel.updateTime()
    numTimeSteps = numTimeSteps + 1
    writer = exporters.VTKWriterA(geomFields,meshes,"flow-" + str(numTimeSteps) + ".vtk",
                              "FlowField",False,0)
    writer.init()
    writer.writeScalarField(flowFields.pressure,"Pressure")
    writer.writeVectorField(flowFields.velocity,"Velocity")
    writer.finish()

t1 = time.time()
if outfile != '/dev/stdout':
    print '\nsolution time = %f' % (t1-t0)

if (atype=='tangent'):
    writer = exporters.FluentDataExporterA(reader,fileBase+"-prism-tangent.dat",False,1)
    writer.init()
    writer.writeScalarField(flowFields.pressure,1)
    writer.writeVectorField(flowFields.velocity,111)
    writer.writeScalarField(flowFields.massFlux,18)
    writer.finish()

    
