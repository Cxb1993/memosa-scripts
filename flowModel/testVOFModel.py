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

rhoAir = 1.225
rhoH20 = 1000.0
muAir = 1.7894E-5
muH20 = 1.002E-3

geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

if atype == 'tangent':
    metricsCalculator.setTangentCoords(0,7,1)

flowFields =  models.FlowFields('flow')
fmodel = models.FlowModelA(geomFields,flowFields,meshes)

thermalFields = models.ThermalFields('temperature')
tmodel = models.ThermalModelA(geomFields,thermalFields,meshes)

reader.importFlowBCs(fmodel,meshes)
bcMap = fmodel.getBCMap()
bcMapTherm = tmodel.getBCMap()
bcMapTherm[7].bcType = 'Symmetry'
#bc = bcMap[7]
#bc.setVar("specifiedXVelocity",0.0)

for i in [5,6,8,9,10]:
    bc = bcMap[i]
    bc.bcType = 'Symmetry'
    #bc.bcType = 'NoSlipWall'
    bcMapTherm[i].bcType = 'Symmetry'
fmodel.printBCs()
tmodel.printBCs()

vcmap = tmodel.getVCMap()
vcZone = vcmap[meshes[0].getID()]
vcZone['thermalConductivity'] = 1.0E-12
vcZone['density'] = rhoAir
vcZone['specificHeat'] = 1.0

momSolver = fvmbaseExt.AMG()
momSolver.relativeTolerance = 1e-15
momSolver.absoluteTolerance = 1.e-15
momSolver.nMaxIterations = 200
momSolver.maxCoarseLevels=20
momSolver.verbosity=0

contSolver = fvmbaseExt.AMG()
#pc = fvmbaseExt.AMG()
#pc.verbosity=0
#contSolver = fvmbaseExt.BCGStab()
#contSolver.preconditioner = pc
contSolver.relativeTolerance = 1e-15
contSolver.absoluteTolerance = 1.e-15
contSolver.nMaxIterations = 200
contSolver.verbosity=0
contSolver.maxCoarseLevels=20

pc = fvmbaseExt.AMG()
pc.verbosity = 0
tempSolver = fvmbaseExt.BCGStab()
tempSolver.preconditioner = pc
tempSolver.relativeTolerance = 1.e-15
tempSolver.absoluteTolerance = 1.e-15
tempSolver.nMaxIteractions = 200
tempSolver.verbosity = 0

foptions = fmodel.getOptions()
toptions = tmodel.getOptions()

foptions.momentumLinearSolver = momSolver
foptions.pressureLinearSolver = contSolver
toptions.linearSolver = tempSolver

foptions.momentumTolerance=1e-3
foptions.continuityTolerance=1e-3
foptions.setVar("momentumURF",0.7)
foptions.setVar("pressureURF",0.3)
foptions.printNormalizedResiduals = True
foptions.transient = True
foptions.setVar("timeStep",0.01)

toptions.relativeTolerance=1e-6
toptions.absoluteTolerance=1e-6
toptions.setVar("initialTemperature",0.0)
toptions.enthalpyModel = False
toptions.polynomialCp = False
toptions.transient = True
toptions.setVar("timeStep",0.01)

"""
if atype=='tangent':
    vcMap = fmodel.getVCMap()
    for i,vc in vcMap.iteritems():
        print vc.getVar('viscosity')
        vc.setVar('viscosity',(1.7894e-5,1))
"""
numTimeSteps = 0
fmodel.init()
tmodel.init()

cells = meshes[0].getCells()
faces = meshes[0].getFaces()
fg = meshes[0].getInteriorFaceGroup()
intFaces = fg.site
faceCells = meshes[0].getAllFaceCells()
intFaceCells = meshes[0].getFaceCells(intFaces)
intFaceArea = geomFields.area[intFaces].asNumPyArray()
CellCoord = geomFields.coordinate[cells].asNumPyArray()

source = flowFields.source[cells].asNumPyArray()
rho = flowFields.density[cells].asNumPyArray()
rhoPrev = flowFields.densityN1[cells].asNumPyArray()
mu = flowFields.viscosity[cells].asNumPyArray()
cellPressure = flowFields.pressure[cells].asNumPyArray()
facePressure = flowFields.pressure[faces].asNumPyArray()
massFlux = flowFields.massFlux[faces].asNumPyArray()
intMassFlux = flowFields.massFlux[intFaces].asNumPyArray()
V = flowFields.velocity[cells].asNumPyArray()

mfrac = thermalFields.temperature[cells].asNumPyArray()
mfracPrev = thermalFields.temperatureN1[cells].asNumPyArray()
cp = thermalFields.specificHeat[cells].asNumPyArray()
cpPrev = thermalFields.specificHeatN1[cells].asNumPyArray()
mfracFlux = thermalFields.convectionFlux[faces].asNumPyArray()

grav = np.array([0.0,-9.81,0.0])

source[:] = 0.0

for c in range(cells.getSelfCount()):
    x = CellCoord[c,0]
    if x <= 0.03:
        mfrac[c] = 1.0
        mfracPrev[c] = 1.0
    rho[c] = mfrac[c]*rhoH20 + (1.0-mfrac[c])*rhoAir
    rhoPrev[c] = rho[c]
    cp[c] = rho[c]
    cpPrev[c] = rho[c]
    source[c] = rho[c]*grav
    mu[c] = mfrac[c]*muH20 + (1.0-mfrac[c])*muAir



pdb.set_trace()
writer = exporters.VTKWriterA(geomFields,meshes,"flow-" + str(numTimeSteps) + ".vtk",
                              "FlowField",False,0)
writer.init()
writer.writeScalarField(flowFields.pressure,"Pressure")
writer.writeScalarField(thermalFields.temperature, "mfrac")
writer.writeScalarField(flowFields.density,"Density")
writer.writeScalarField(flowFields.viscosity,"Viscosity")
writer.writeVectorField(flowFields.velocity,"Velocity")
writer.finish()

while (numTimeSteps < 100):
    for j in range(300):
        flowConverged = fmodel.advance(1)
        mfracFlux[:] = massFlux[:]
        mfracConverged = tmodel.advance(1)
        #Update density and source fields
        for c in range(cells.getSelfCount()):
            vof = mfrac[c]*rho[c]/rhoH20                                                                                   
            rho[c] = vof*rhoH20 + (1.0-vof)*rhoAir
            cp[c] = rho[c]                                                                                                    
            mu[c] = vof*muH20 + (1.0-vof)*muAir                                                                               
            source[c] = grav*rho[c]
            
        #Update mass flux at interior faces
        for f in range(intFaces.getCount()):
            c0 = intFaceCells(f,0)
            c1 = intFaceCells(f,1)
            intMassFlux[f] = 0.5*(np.dot(V[c0,:],intFaceArea[f,:])*rho[c0] +
                                  np.dot(V[c1,:],intFaceArea[f,:])*rho[c1])

        if(flowConverged and mfracConverged):
            break
    else:
        raise NameError('Model did not converge in given number of iterations')
    fmodel.updateTime()
    tmodel.updateTime()
    numTimeSteps = numTimeSteps + 1
    writer = exporters.VTKWriterA(geomFields,meshes,"flow-" + str(numTimeSteps) + ".vtk",
                              "FlowField",False,0)
    writer.init()
    writer.writeScalarField(flowFields.pressure,"Pressure")
    writer.writeScalarField(thermalFields.temperature, "mfrac")
    writer.writeScalarField(flowFields.density,"Density")
    writer.writeScalarField(flowFields.viscosity,"Viscosity")
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

    
