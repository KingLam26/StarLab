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

# =============== L12542-step2-part1 =============== #
# set up parameters
LensOneFullHeight = 12.7
LensOneRadius = 8.7
LensOneEdgeThickness = 1.5
LensOneCentreThickness = functions.calculatePCLensCT(LensOneFullHeight, LensOneRadius, LensOneEdgeThickness)

# define function
def step2part1(SourcesLensDistance):
    global beam1, beam2, beam3
    global LensOnePos1, LensOnePos2

    sld = SourcesLensDistance[0]

    beam1 = beamClass.Beam([0, 0.3], np.radians(12), 'beam1')
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
print(f'\nLensOnePos1: {LensOnePos1}')
print(f'LensOnePos2: {LensOnePos2}')
print(f'SourcesLensDistance: {SourcesLensDistance}')
print(f'SourcesLensDistance: {SourcesLensDistance.x[0]}')


# =============== L12542-step2-part2 =============== #
# set up parameters
LensTwoFullHeight = 25.4
LensTwoRadius = 216.9
LensTwoEdgeThickness = 2
LensTwoCentreThickness = functions.calculatePCLensCT(LensTwoFullHeight, LensTwoRadius, LensTwoEdgeThickness)

LensOneTwoDistance = 30
LensTwoEntryEquation = [SourcesLensDistance.x[0] + LensOneTwoDistance + LensTwoRadius - LensTwoCentreThickness, 0, LensTwoRadius]
LensTwoPos1 = beam3.propagateBeam('circle', LensTwoEntryEquation, closestFlag=True)
refractedBeam = beam3.refractBeam('circle', LensTwoEntryEquation, closestFlag=True, RI1 = 1, RI2 = LensRI)
beam4 = beamClass.Beam(LensTwoPos1, functions.lineAngle(refractedBeam), 'beam4')

LensTwoPos2 = beam4.propagateBeam('line', [1,0,-(SourcesLensDistance.x[0] + LensOneTwoDistance)], closestFlag=True)
refractedBeam = beam4.refractBeam('line', [1,0,-(SourcesLensDistance.x[0] + LensOneTwoDistance)], closestFlag=True,RI1 = LensRI, RI2 = 1)
beam5 = beamClass.Beam(LensTwoPos2, functions.lineAngle(refractedBeam), 'beam5')

beam3.printBeam()
beam4.printBeam()
print(f'\nLensTwoPos1: {LensTwoPos1}')
print(f'LensTwoPos2: {LensTwoPos2}')
print(f'LensOneTwoDistance: {LensOneTwoDistance}')


# =============== L12542-step2-part3 =============== #
# set up parameters: prism related
PrismLength = 25.4
PrismTilt = np.radians(0)

LensTwoPrismDistance = 50       # variable

def step2part3(LensTwoPrismDistance):
    global beam5, beam6, beam7, beam8, beam9, beam10

    LensTwoPrismDistance = LensTwoPrismDistance[0]
    
    PrismPosition = [SourcesLensDistance.x[0] + LensOneTwoDistance + LensTwoPrismDistance, 0]
    PrismEntryCorner1 = [PrismPosition[0] - np.sin(PrismTilt)*PrismLength/2, PrismPosition[1] + np.cos(PrismTilt)*PrismLength/2]
    PrismEntryCorner2 = [PrismPosition[0] + np.sin(PrismTilt)*PrismLength/2, PrismPosition[1] - np.cos(PrismTilt)*PrismLength/2]
    PrismExitCorner1 = PrismEntryCorner2
    PrismExitCorner2 = [PrismExitCorner1[0] + np.cos(PrismTilt)*PrismLength, PrismExitCorner1[1] - np.sin(PrismTilt)*PrismLength]

    PrismEntryLineEqn = functions.lineEquationPoints(PrismEntryCorner1, PrismEntryCorner2)
    PrismExitLineEqn = functions.lineEquationPoints(PrismExitCorner1, PrismExitCorner2)
    PrismHypotenusEqn = functions.lineEquationPoints(PrismEntryCorner1, PrismExitCorner2)

    PrismPos1 = beam5.propagateBeam('line', PrismEntryLineEqn, closestFlag=True)
    refractedBeam = beam5.refractBeam('line', PrismEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = PrismRI)
    beam6 = beamClass.Beam(PrismPos1, functions.lineAngle(refractedBeam), 'beam6')
    PrismPos2 = beam6.propagateBeam('line', PrismHypotenusEqn, closestFlag=True)

    reflectedBeam = beam6.reflectBeam('line', PrismHypotenusEqn, closestFlag=True)
    beam7 = beamClass.Beam(PrismPos2, functions.lineAngle(reflectedBeam), 'beam7')
    PrismPos3 = beam7.propagateBeam('line', PrismExitLineEqn, closestFlag=True)

    refractedBeam = beam7.refractBeam('line', PrismExitLineEqn, closestFlag=True, RI1 = PrismRI, RI2 = 1)
    beam8 = beamClass.Beam(PrismPos3, functions.lineAngle(refractedBeam), 'beam8')

    # set up parameters: prism to sample
    prismSampleDistance = 400
    topWindowSampleDistance = 165.3
    windowThickmess = 2.997
    prismTopWindowDistance = prismSampleDistance - topWindowSampleDistance
    bottomWindowSampleDistance = prismTopWindowDistance - windowThickmess

    windowPosition = [PrismPosition[0] + PrismLength/2, PrismPosition[1] -PrismLength/2 - prismTopWindowDistance]
    windowEntryLineEqn = [0, 1, -windowPosition[1]]
    windowExitLineEqn = [0, 1, -(windowPosition[1] - windowThickmess)]

    targetHoriztonalLineEqn = [0, 1, -(windowPosition[1]-windowThickmess-bottomWindowSampleDistance)]
    targetPosition = [windowPosition[0], windowPosition[1]-windowThickmess-bottomWindowSampleDistance]

    windowPos1 = beam8.propagateBeam('line', windowEntryLineEqn, closestFlag=True )
    refractedBeam = beam8.refractBeam('line', windowEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = WindowsRI)
    beam9 = beamClass.Beam(windowPos1, functions.lineAngle(refractedBeam), 'beam9')
    windowPos2 = beam9.propagateBeam('line', windowExitLineEqn, closestFlag=True)
    refractedBeam = beam9.refractBeam('line', windowExitLineEqn, closestFlag=True, RI1 = WindowsRI, RI2 = 1)
    beam10 = beamClass.Beam(windowPos2, functions.lineAngle(refractedBeam), 'beam10')
    samplePosition = beam10.propagateBeam('line', targetHoriztonalLineEqn, closestFlag=True)

    targetSampleDistance = functions.pointPointDistance(samplePosition, targetPosition)

    return targetSampleDistance

LensTwoPrismDistance = minimize(step2part3, 2, tol=1e-10, method = 'Nelder-Mead')

print(f'LensTwoPrismDistance: {LensTwoPrismDistance}')
print(f'LensTwoPrismDistance: {LensTwoPrismDistance.x[0]}')


beam5.printBeam()
beam6.printBeam()
beam7.printBeam()
beam8.printBeam()
beam9.printBeam()
beam10.printBeam()