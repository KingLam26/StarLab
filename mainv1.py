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
    global beamA1, beamA2, beamA3
    global LensOnePosA1, LensOnePosA2

    sld = SourcesLensDistance[0]

    beamA1 = beamClass.Beam([0, 0.3], np.radians(12), 'beamA1')
    LensOnePosA1 = beamA1.propagateBeam('line', [1,0,-sld], closestFlag=True)
    beamA1.refractBeam('line', [1,0,-sld], closestFlag=True, RI1 = 1, RI2 = LensRI)

    beamA2 = beamClass.Beam(LensOnePosA1, beamA1.ExitAngle, 'beamA2')
    LensOnePosA2 = beamA2.propagateBeam('circle', [sld - (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True)
    refracted = beamA2.refractBeam('circle', [sld - (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True, RI1 = LensRI, RI2 = 1)

    beamA3 = beamClass.Beam(LensOnePosA2, functions.lineAngle(refracted), 'beamA3')
    return np.abs(beamA3.gradient-0)

SourcesLensDistance = minimize(step2part1, 2, tol=1e-10, method = 'Nelder-Mead')

# print information about the beams
print(f'\nLensOnePosA1: {LensOnePosA1}')
print(f'LensOnePosA2: {LensOnePosA2}')
print(f'SourcesLensDistance: {SourcesLensDistance}')
print(f'SourcesLensDistance: {SourcesLensDistance.x[0]}')


# =============== L12542-step2-part2 =============== #
# set up parameters
LensTwoFullHeight = 25.4
LensTwoRadius = 216.9
LensTwoEdgeThickness = 2
LensTwoCentreThickness = functions.calculatePCLensCT(LensTwoFullHeight, LensTwoRadius, LensTwoEdgeThickness)

LensOneTwoDistance = 50
LensTwoEntryEquation = [SourcesLensDistance.x[0] + LensOneTwoDistance + LensTwoRadius - LensTwoCentreThickness, 0, LensTwoRadius]
LensTwoPosA1 = beamA3.propagateBeam('circle', LensTwoEntryEquation, closestFlag=True)
refractedBeam = beamA3.refractBeam('circle', LensTwoEntryEquation, closestFlag=True, RI1 = 1, RI2 = LensRI)
beamA4 = beamClass.Beam(LensTwoPosA1, functions.lineAngle(refractedBeam), 'beamA4')

LensTwoPosA2 = beamA4.propagateBeam('line', [1,0,-(SourcesLensDistance.x[0] + LensOneTwoDistance)], closestFlag=True)
refractedBeam = beamA4.refractBeam('line', [1,0,-(SourcesLensDistance.x[0] + LensOneTwoDistance)], closestFlag=True,RI1 = LensRI, RI2 = 1)
beamA5 = beamClass.Beam(LensTwoPosA2, functions.lineAngle(refractedBeam), 'beamA5')

print(f'\nLensTwoPosA1: {LensTwoPosA1}')
print(f'LensTwoPosA2: {LensTwoPosA2}')
print(f'LensOneTwoDistance: {LensOneTwoDistance}')


# =============== L12542-step2-part3 =============== #
# set up parameters: prism related
PrismLength = 25.4
PrismTilt = np.radians(0)

LensTwoPrismDistance = 50       # variable

def step2part3(LensTwoPrismDistance):
    global beamA5, beamA6, beamA7, beamA8, beamA9, beamA10
    global targetSampleDistance, samplePosition, targetPosition

    LensTwoPrismDistance = LensTwoPrismDistance[0]
    
    PrismPosition = [SourcesLensDistance.x[0] + LensOneTwoDistance + LensTwoPrismDistance, 0]
    PrismEntryCorner1 = [PrismPosition[0] - np.sin(PrismTilt)*PrismLength/2, PrismPosition[1] + np.cos(PrismTilt)*PrismLength/2]
    PrismEntryCorner2 = [PrismPosition[0] + np.sin(PrismTilt)*PrismLength/2, PrismPosition[1] - np.cos(PrismTilt)*PrismLength/2]
    PrismExitCorner1 = PrismEntryCorner2
    PrismExitCorner2 = [PrismExitCorner1[0] + np.cos(PrismTilt)*PrismLength, PrismExitCorner1[1] - np.sin(PrismTilt)*PrismLength]

    PrismEntryLineEqn = functions.lineEquationPoints(PrismEntryCorner1, PrismEntryCorner2)
    PrismExitLineEqn = functions.lineEquationPoints(PrismExitCorner1, PrismExitCorner2)
    PrismHypotenusEqn = functions.lineEquationPoints(PrismEntryCorner1, PrismExitCorner2)

    PrismPos1 = beamA5.propagateBeam('line', PrismEntryLineEqn, closestFlag=True)
    refractedBeam = beamA5.refractBeam('line', PrismEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = PrismRI)
    beamA6 = beamClass.Beam(PrismPos1, functions.lineAngle(refractedBeam), 'beamA6')
    PrismPos2 = beamA6.propagateBeam('line', PrismHypotenusEqn, closestFlag=True)

    reflectedBeam = beamA6.reflectBeam('line', PrismHypotenusEqn, closestFlag=True)
    beamA7 = beamClass.Beam(PrismPos2, functions.lineAngle(reflectedBeam), 'beamA7')
    PrismPos3 = beamA7.propagateBeam('line', PrismExitLineEqn, closestFlag=True)

    refractedBeam = beamA7.refractBeam('line', PrismExitLineEqn, closestFlag=True, RI1 = PrismRI, RI2 = 1)
    beamA8 = beamClass.Beam(PrismPos3, functions.lineAngle(refractedBeam), 'beamA8')

    # set up parameters: prism to sample
    prismSampleDistance = 430
    topWindowSampleDistance = 165.3
    windowThickmess = 2.997
    prismTopWindowDistance = prismSampleDistance - topWindowSampleDistance
    bottomWindowSampleDistance = prismSampleDistance - windowThickmess - prismTopWindowDistance

    windowPosition = [PrismPosition[0] + PrismLength/2, PrismPosition[1] -PrismLength/2 - prismTopWindowDistance]
    windowEntryLineEqn = [0, 1, -windowPosition[1]]
    windowExitLineEqn = [0, 1, -(windowPosition[1] - windowThickmess)]

    targetHoriztonalLineEqn = [0, 1, -(windowPosition[1]-windowThickmess-bottomWindowSampleDistance)]
    targetPosition = [windowPosition[0], windowPosition[1]-windowThickmess-bottomWindowSampleDistance]

    windowPosA1 = beamA8.propagateBeam('line', windowEntryLineEqn, closestFlag=True )
    refractedBeam = beamA8.refractBeam('line', windowEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = WindowsRI)
    beamA9 = beamClass.Beam(windowPosA1, functions.lineAngle(refractedBeam), 'beamA9')
    windowPosA2 = beamA9.propagateBeam('line', windowExitLineEqn, closestFlag=True)
    refractedBeam = beamA9.refractBeam('line', windowExitLineEqn, closestFlag=True, RI1 = WindowsRI, RI2 = 1)
    beamA10 = beamClass.Beam(windowPosA2, functions.lineAngle(refractedBeam), 'beamA10')
    samplePosition = beamA10.propagateBeam('line', targetHoriztonalLineEqn, closestFlag=True)

    targetSampleDistance = functions.pointPointDistance(samplePosition, targetPosition)

    return targetSampleDistance

LensTwoPrismDistance = minimize(step2part3, 2, tol=1e-10, method = 'Nelder-Mead')





beamA1.printBeam()
beamA2.printBeam()
beamA3.printBeam()
beamA4.printBeam()
beamA5.printBeam()
beamA6.printBeam()
beamA7.printBeam()
beamA8.printBeam()
beamA9.printBeam()
beamA10.printBeam()

print(f'LensTwoPrismDistance: {LensTwoPrismDistance}')
print(f'LensTwoPrismDistance: {LensTwoPrismDistance.x[0]}')

print(f'samplePosition: {samplePosition}')
print(f'targetPosition: {samplePosition}')
