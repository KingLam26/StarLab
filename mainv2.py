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
# propagate from source to L1
# set up parameters
LensOneFullHeight = 25.4
LensOneRadius = 21.7
LensOneEdgeThickness = 2
LensOneCentreThickness = functions.calculatePCLensCT(LensOneFullHeight, LensOneRadius, LensOneEdgeThickness, type = 'convex')

sld = 50

beamA1 = beamClass.Beam([0, 0], np.radians(12), 'beamA1')
LensOnePosA1 = beamA1.propagateBeam('line', [1,0,-sld], closestFlag=True)
beamA1.refractBeam('line', [1,0,-sld], closestFlag=True, RI1 = 1, RI2 = LensRI)

beamA2 = beamClass.Beam(LensOnePosA1, beamA1.ExitAngle, 'beamA2')
LensOnePosA2 = beamA2.propagateBeam('circle', [sld - (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True)
refracted = beamA2.refractBeam('circle', [sld - (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True, RI1 = LensRI, RI2 = 1)

beamA3 = beamClass.Beam(LensOnePosA2, functions.lineAngle(refracted), 'beamA3')


# =============== L12542-step2-part2 =============== #
# propagate from L1 to L2
# set up parameters
LensTwoFullHeight = 25.4
LensTwoRadius = 32.5
LensTwoEdgeThickness = 5.1
LensTwoCentreThickness = functions.calculatePCLensCT(LensTwoFullHeight, LensTwoRadius, LensTwoEdgeThickness, type ='concave')

def step2part2(variable):
    
    global beamA4, beamA5

    LensOneTwoDistance = variable[0]     # variable

    LensTwoEntryEquation = [1, 0, -(sld + LensOneTwoDistance)]
    LensTwoPosA1 = beamA3.propagateBeam('line', LensTwoEntryEquation, closestFlag=True)
    refractedBeam = beamA3.refractBeam('line', LensTwoEntryEquation, closestFlag=True, RI1 = 1, RI2 = LensRI)
    beamA4 = beamClass.Beam(LensTwoPosA1, functions.lineAngle(refractedBeam), 'beamA4')

    LensTwoExitEquation = [sld + LensOneTwoDistance + LensTwoRadius + LensTwoCentreThickness, 0, LensTwoRadius]

    LensTwoPosA2 = beamA4.propagateBeam('circle', LensTwoExitEquation, closestFlag=True)
    refractedBeam = beamA4.refractBeam('circle', LensTwoExitEquation, closestFlag=True,RI1 = LensRI, RI2 = 1)
    beamA5 = beamClass.Beam(LensTwoPosA2, functions.lineAngle(refractedBeam), 'beamA5')

    return np.abs(beamA5.gradient)

LensOneTwoDistance = minimize(step2part2, 2, tol=1e-10, method = 'Nelder-Mead')

print(f'LensOneTwoDistance: {LensOneTwoDistance}')


# =============== L12542-step2-part3 =============== #
# propagate from L2 to L3

LensTwoThreeDistance = 50

LensThreeFullHeight = 25.4
LensThreeRadius = 216.9
LensThreeEdgeThickness = 2
LensThreeCentreThickness = functions.calculatePCLensCT(LensThreeFullHeight, LensThreeRadius, LensThreeEdgeThickness, 'convex')

LensThreeEntryEquation = [(sld + LensOneTwoDistance.x[0] + LensTwoThreeDistance + LensThreeRadius - LensThreeCentreThickness), 0, LensThreeRadius]

LensThreePosA1 = beamA5.propagateBeam('circle', LensThreeEntryEquation, closestFlag=True)
refractedBeam = beamA5.refractBeam('circle', LensThreeEntryEquation, closestFlag=True, RI1 = 1, RI2 = LensRI)
beamA6 = beamClass.Beam(LensThreePosA1, functions.lineAngle(refractedBeam), 'beamA6')

LensThreeExitEquation = [1, 0, -(sld + LensOneTwoDistance.x[0] + LensTwoThreeDistance)]

LensThreePosA2 = beamA6.propagateBeam('line', LensThreeExitEquation, closestFlag=True)
refractedBeam = beamA6.refractBeam('line', LensThreeExitEquation, closestFlag=True,RI1 = LensRI, RI2 = 1)
beamA7 = beamClass.Beam(LensThreePosA2, functions.lineAngle(refractedBeam), 'beamA7')


# =============== L12542-step2-part3 =============== #
# propagate from L3 to prism
LensThreePrismDistance = 51.3       # variable

PrismLength = 25.4
PrismTilt = np.radians(0)


def step2part4(variable):
    global beamA8, beamA9, beamA10, beamA11, beamA12
    global targetSampleDistance, samplePosition, targetPosition

    LensThreePrismDistance = variable[0]

    PrismPosition = [sld + LensOneTwoDistance.x[0] + LensTwoThreeDistance + LensThreePrismDistance, 0]
    PrismEntryCorner1 = [PrismPosition[0] - np.sin(PrismTilt)*PrismLength/2, PrismPosition[1] + np.cos(PrismTilt)*PrismLength/2]
    PrismEntryCorner2 = [PrismPosition[0] + np.sin(PrismTilt)*PrismLength/2, PrismPosition[1] - np.cos(PrismTilt)*PrismLength/2]
    PrismExitCorner1 = PrismEntryCorner2
    PrismExitCorner2 = [PrismExitCorner1[0] + np.cos(PrismTilt)*PrismLength, PrismExitCorner1[1] - np.sin(PrismTilt)*PrismLength]

    PrismEntryLineEqn = functions.lineEquationPoints(PrismEntryCorner1, PrismEntryCorner2)
    PrismExitLineEqn = functions.lineEquationPoints(PrismExitCorner1, PrismExitCorner2)
    PrismHypotenusEqn = functions.lineEquationPoints(PrismEntryCorner1, PrismExitCorner2)

    PrismPos1 = beamA7.propagateBeam('line', PrismEntryLineEqn, closestFlag=True)
    refractedBeam = beamA7.refractBeam('line', PrismEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = PrismRI)
    beamA8 = beamClass.Beam(PrismPos1, functions.lineAngle(refractedBeam), 'beamA8')
    PrismPos2 = beamA8.propagateBeam('line', PrismHypotenusEqn, closestFlag=True)

    reflectedBeam = beamA8.reflectBeam('line', PrismHypotenusEqn, closestFlag=True)
    beamA9 = beamClass.Beam(PrismPos2, functions.lineAngle(reflectedBeam), 'beamA9')
    PrismPos3 = beamA9.propagateBeam('line', PrismExitLineEqn, closestFlag=True)
    refractedBeam = beamA9.refractBeam('line', PrismExitLineEqn, closestFlag=True, RI1 = PrismRI, RI2 = 1)
    beamA10 = beamClass.Beam(PrismPos3, functions.lineAngle(refractedBeam), 'beamA10')

    # propagate from prism to sample
    prismSampleDistance = 400
    topWindowSampleDistance = 165.3
    windowThickmess = 2.997
    prismTopWindowDistance = prismSampleDistance - topWindowSampleDistance
    bottomWindowSampleDistance = prismSampleDistance - windowThickmess - prismTopWindowDistance

    windowPosition = [PrismPosition[0] + PrismLength/2, PrismPosition[1] -PrismLength/2 - prismTopWindowDistance]
    windowEntryLineEqn = [0, 1, -windowPosition[1]]
    windowExitLineEqn = [0, 1, -(windowPosition[1] - windowThickmess)]

    targetHoriztonalLineEqn = [0, 1, -(windowPosition[1]-windowThickmess-bottomWindowSampleDistance)]
    targetPosition = [windowPosition[0], windowPosition[1]-windowThickmess-bottomWindowSampleDistance]

    windowPosA1 = beamA10.propagateBeam('line', windowEntryLineEqn, closestFlag=True )
    refractedBeam = beamA10.refractBeam('line', windowEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = WindowsRI)
    beamA11 = beamClass.Beam(windowPosA1, functions.lineAngle(refractedBeam), 'beamA11')
    windowPosA2 = beamA11.propagateBeam('line', windowExitLineEqn, closestFlag=True)
    refractedBeam = beamA11.refractBeam('line', windowExitLineEqn, closestFlag=True, RI1 = WindowsRI, RI2 = 1)
    beamA12 = beamClass.Beam(windowPosA2, functions.lineAngle(refractedBeam), 'beamA12')
    samplePosition = beamA12.propagateBeam('line', targetHoriztonalLineEqn, closestFlag=True)

    targetSampleDistance = functions.pointPointDistance(samplePosition, targetPosition)

    return targetSampleDistance

LensThreePrismDistance = minimize(step2part4, 2, tol=1e-10, method = 'Nelder-Mead')


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
beamA11.printBeam()
beamA12.printBeam()


print(f'LensTwoPrismDistance: {LensThreePrismDistance}')

print(f'samplePosition: {samplePosition}')
print(f'targetPosition: {samplePosition}')
