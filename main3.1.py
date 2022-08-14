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

# =============== EQ-step2-part1 =============== #
# set up parameters
sysTwoOrigin = []
LensOneFullHeight = 12.7
LensOneRadius = 8.7
LensOneEdgeThickness = 1.5
LensOneCentreThickness = functions.calculatePCLensCT(LensOneFullHeight, LensOneRadius, LensOneEdgeThickness)
clearApertureRadius = 0.9*LensOneFullHeight/2

# define function
def step2part1(SourcesLensDistance):
    global beam1, beam2, beam3

    global LensOnePos1

    sld = SourcesLensDistance[0]

    beam1 = beamClass.Beam([0, 0], np.radians(150), 'beam1')
    LensOnePos1 = beam1.propagateBeam('line', [1,0,sld], closestFlag=True)
    beam1.refractBeam('line', [1,0,sld], closestFlag=True, RI1 = 1, RI2 = LensRI)

    beam2 = beamClass.Beam(LensOnePos1, np.pi - beam1.ExitAngle, 'beam2')
    LensOnePos2 = beam2.propagateBeam('circle', [-sld + (LensOneRadius - LensOneCentreThickness), 0, LensOneRadius], closestFlag=True)
    refracted = beam2.refractBeam('circle', [-sld + (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True, RI1 = LensRI, RI2 = 1)

    beam3 = beamClass.Beam(LensOnePos2, functions.lineAngle(refracted), 'beam3')

    return np.abs(LensOnePos2[1] - clearApertureRadius)

SourcesLensDistance = minimize(step2part1, 2, tol=1e-10, method = 'Nelder-Mead')

# print information about the beams
beam1.printBeam()
beam2.printBeam()
print(f'SourcesLensDistance: {SourcesLensDistance.x[0]}')


# =============== EQ-step2-part2 =============== #
# set up parameters
LensTwoFullHeight = 25.4
LensTwoRadius = 216.9
LensTwoEdgeThickness = 2
LensTwoCentreThickness = functions.calculatePCLensCT(LensTwoFullHeight, LensTwoRadius, LensTwoEdgeThickness)

LensOneTwoDistance = 10                 # variable 1
LensTwoEntryEquation = [-SourcesLensDistance.x[0] - LensOneTwoDistance - LensTwoRadius + LensTwoCentreThickness, 0, LensTwoRadius]

LensTwoPos1 = beam3.propagateBeam('circle', LensTwoEntryEquation, closestFlag=True)

refractedBeam = beam3.refractBeam('circle', LensTwoEntryEquation, closestFlag=True, RI1 = 1, RI2 = LensRI)
beam4 = beamClass.Beam(LensTwoPos1, functions.lineAngle(refractedBeam), 'beam4')

LensTwoPos2 = beam4.propagateBeam('line', [1,0,(SourcesLensDistance.x[0] + LensOneTwoDistance)], closestFlag=True)
refractedBeam = beam4.refractBeam('line', [1,0,(SourcesLensDistance.x[0] + LensOneTwoDistance)], closestFlag=True,RI1 = LensRI, RI2 = 1)
beam5 = beamClass.Beam(LensTwoPos2, functions.lineAngle(refractedBeam), 'beam5')

beam3.printBeam()
beam4.printBeam()
print(f'refractedBeam:{refractedBeam}')
print(beam5.gradient)