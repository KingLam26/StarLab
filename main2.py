# =============== import libraries =============== #
from fileinput import close
from math import radians
from re import T
import numpy as np
from scipy.optimize import minimize
import functions, beamClass

# =============== general parameters =============== #
WindowsRI = 1.459
LensRI = 1.434
PrismRI = 1.434

# =============== Horiba-step2-part1 =============== #
# set up parameters
sysTwoOrigin = []
LensOneFullHeight = 12.7
LensOneRadius = 8.7
LensOneEdgeThickness = 1.5
LensOneCentreThickness = functions.calculatePCLensCT(LensOneFullHeight, LensOneRadius, LensOneEdgeThickness)

# define function
def step2part1(SourcesLensDistance):
    global beam1, beam2, beam3
    global LensOnePos1, LensOnePos2

    sld = SourcesLensDistance[0]

    beam1 = beamClass.Beam([0, 0], np.radians(14), 'beam1')
    LensOnePos1 = beam1.propagateBeam('line', [1,0,-sld], closestFlag=True)
    beam1.refractBeam('line', [1,0,-sld], closestFlag=True, RI1 = 1, RI2 = LensRI)

    beam2 = beamClass.Beam(LensOnePos1, beam1.ExitAngle, 'beam2')
    LensOnePos2 = beam2.propagateBeam('circle', [sld - (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True)
    refracted = beam2.refractBeam('circle', [sld - (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True, RI1 = LensRI, RI2 = 1)

    beam3 = beamClass.Beam(LensOnePos2, functions.lineAngle(refracted), 'beam3')
    return np.abs(beam3.gradient-0)

SourcesLensDistance = minimize(step2part1, 2, tol=1e-10, method = 'Nelder-Mead')

# print information about the beams
beam1.printBeam()
beam2.printBeam()
beam3.printBeam()
print(f'\nLensOnePos1: {LensOnePos1}')
print(f'LensOnePos2: {LensOnePos2}')
print(f'SourcesLensDistance: {SourcesLensDistance.x[0]}')