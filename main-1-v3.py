"""
For L12542 / Deuterium lamp
"""

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
# MirrorRI = 1.434
MirrorRI = 1

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
print(f'LensTwoCentreThickness: {LensTwoCentreThickness}')

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
# propagate from L3 to mirror
LensThreeMirrorDistance = 51.3       # variable

MirrorHypotenuse = 45
MirrorLength = np.sqrt(MirrorHypotenuse**2 / 2)
print(f'MirrorLength: {MirrorLength}')
MirrorTilt = np.radians(0)


def step2part4(variable):
    global beamA8, beamA9, beamA10, beamA11, beamA12
    global targetSampleDistance, samplePosition, targetPosition

    LensThreeMirrorDistance = variable[0]

    MirrorPosition = [sld + LensOneTwoDistance.x[0] + LensTwoThreeDistance + LensThreeMirrorDistance, 0]
    MirrorEntryCorner1 = [MirrorPosition[0] - np.sin(MirrorTilt)*MirrorLength/2, MirrorPosition[1] + np.cos(MirrorTilt)*MirrorLength/2]
    MirrorEntryCorner2 = [MirrorPosition[0] + np.sin(MirrorTilt)*MirrorLength/2, MirrorPosition[1] - np.cos(MirrorTilt)*MirrorLength/2]
    MirrorExitCorner1 = MirrorEntryCorner2
    MirrorExitCorner2 = [MirrorExitCorner1[0] + np.cos(MirrorTilt)*MirrorLength, MirrorExitCorner1[1] - np.sin(MirrorTilt)*MirrorLength]

    MirrorEntryLineEqn = functions.lineEquationPoints(MirrorEntryCorner1, MirrorEntryCorner2)
    MirrorExitLineEqn = functions.lineEquationPoints(MirrorExitCorner1, MirrorExitCorner2)
    MirrorHypotenusEqn = functions.lineEquationPoints(MirrorEntryCorner1, MirrorExitCorner2)

    MirrorPos1 = beamA7.propagateBeam('line', MirrorEntryLineEqn, closestFlag=True)
    refractedBeam = beamA7.refractBeam('line', MirrorEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = MirrorRI)
    beamA8 = beamClass.Beam(MirrorPos1, functions.lineAngle(refractedBeam), 'beamA8')
    MirrorPos2 = beamA8.propagateBeam('line', MirrorHypotenusEqn, closestFlag=True)

    reflectedBeam = beamA8.reflectBeam('line', MirrorHypotenusEqn, closestFlag=True)
    beamA9 = beamClass.Beam(MirrorPos2, functions.lineAngle(reflectedBeam), 'beamA9')
    MirrorPos3 = beamA9.propagateBeam('line', MirrorExitLineEqn, closestFlag=True)
    refractedBeam = beamA9.refractBeam('line', MirrorExitLineEqn, closestFlag=True, RI1 = MirrorRI, RI2 = 1)
    beamA10 = beamClass.Beam(MirrorPos3, functions.lineAngle(refractedBeam), 'beamA10')

    # propagate from Mirror to sample
    MirrorSampleDistance = 425
    topWindowSampleDistance = 165.3
    windowThickmess = 2.997
    MirrorTopWindowDistance = MirrorSampleDistance - topWindowSampleDistance
    bottomWindowSampleDistance = MirrorSampleDistance - windowThickmess - MirrorTopWindowDistance

    windowPosition = [MirrorPosition[0] + MirrorLength/2, MirrorPosition[1] -MirrorLength/2 - MirrorTopWindowDistance]
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

LensThreeMirrorDistance = minimize(step2part4, 2, tol=1e-10, method = 'Nelder-Mead')


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


print(f'LensTwoMirrorDistance: {LensThreeMirrorDistance}')

print(f'samplePosition: {samplePosition}')
print(f'targetPosition: {samplePosition}')
