"""
For EQ-77
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
LensRI = 1.4585
# MirrorRI = 1.434
MirrorRI = 1

# =============== EQ - step 1 =============== #
""" Propagate light from EQ-77 source to first lens
"""

# set up parameters
sysTwoOrigin = []
LensOneFullHeight = 75
LensOneRadius = 41.4
LensOneEdgeThickness = 3
LensOneCentreThickness = functions.calculatePCLensCT(LensOneFullHeight, LensOneRadius, LensOneEdgeThickness, type = 'convex')
print(f'LensOneCentreThickness: {LensOneCentreThickness}')
clearAperture = 0.9*LensOneFullHeight

# propagate to first lens
sld = ( clearAperture / 2 ) / np.tan(np.radians(30))
print(sld)
beamA1 = beamClass.Beam([0, 0], np.radians(30), 'beamA1')
LensOnePosA1 = beamA1.propagateBeam('line', [1,0,-sld], closestFlag=True)
beamA1.refractBeam('line', [1,0,-sld], closestFlag=True, RI1 = 1, RI2 = LensRI)

beamA2 = beamClass.Beam(LensOnePosA1, beamA1.ExitAngle, 'beamA2')
LensOnePosA2 = beamA2.propagateBeam('circle', [sld - (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True)
refracted = beamA2.refractBeam('circle', [sld - (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True, RI1 = LensRI, RI2 = 1)

beamA3 = beamClass.Beam(LensOnePosA2, functions.lineAngle(refracted), 'beamA3')


# =============== EQ - step 2 =============== #
""" Propagate light from first lens to second lens
"""
# set up parameters
LensTwoFullHeight = 25.4
LensTwoRadius = 17.4
LensTwoEdgeThickness = 7.3
LensTwoCentreThickness = functions.calculatePCLensCT(LensTwoFullHeight, LensTwoRadius, LensTwoEdgeThickness, type = 'concave')
print(f'LensTwoCentreThickness: {LensTwoCentreThickness}')


# define function
def step2(variable):
    global beamA4, beamA5, beamA6, sld, LensTwoEntryEquation, LensTwoExitEquation

    LensOneTwoDistance = variable[0]     # variable
    # print(f"LensOneTwoDistance: {LensOneTwoDistance}")

    LensTwoEntryEquation = [1, 0, -(sld + LensOneTwoDistance)]
    LensTwoPosA1 = beamA3.propagateBeam('line', LensTwoEntryEquation, closestFlag=True)
    refractedBeam = beamA3.refractBeam('line', LensTwoEntryEquation, closestFlag=True, RI1 = 1, RI2 = LensRI)
    beamA4 = beamClass.Beam(LensTwoPosA1, functions.lineAngle(refractedBeam), 'beamA4')

    LensTwoExitEquation = [sld + LensOneTwoDistance + LensTwoRadius + LensTwoCentreThickness, 0, LensTwoRadius]

    LensTwoPosA2 = beamA4.propagateBeam('circle', LensTwoExitEquation, closestFlag=True)
    refractedBeam = beamA4.refractBeam('circle', LensTwoExitEquation, closestFlag=True,RI1 = LensRI, RI2 = 1)
    beamA5 = beamClass.Beam(LensTwoPosA2, functions.lineAngle(refractedBeam), 'beamA5')

    # print(f"beamA5.gradient: {beamA5.gradient}")
    return np.abs(beamA5.gradient)

LensOneTwoDistance = minimize(step2, 200, tol=1e-10, method = 'Nelder-Mead')

# print information about the beams
print(f'LensOneTwoDistance: {LensOneTwoDistance.x[0]}')



# =============== EQ - step 3 =============== #
""" Propagate light from second lens to third lens
"""

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


# =============== EQ-77 step 4 =============== #
""" construct beam B from source to L3
"""

beamB1 = beamClass.Beam([0, 0], np.radians(-30), 'beamB1')
LensOnePosB1 = beamB1.propagateBeam('line', [1,0,-sld], closestFlag=True)
beamB1.refractBeam('line', [1,0,-sld], closestFlag=True, RI1 = 1, RI2 = LensRI)

beamB2 = beamClass.Beam(LensOnePosB1, -beamB1.ExitAngle, 'beamB2')
LensOnePosB2 = beamB2.propagateBeam('circle', [sld - (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True)
refracted = beamB2.refractBeam('circle', [sld - (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True, RI1 = LensRI, RI2 = 1)

beamB3 = beamClass.Beam(LensOnePosB2, functions.lineAngle(refracted), 'beamB3')
LensTwoPosB1 = beamB3.propagateBeam('line', LensTwoEntryEquation, closestFlag=True)
refractedBeam = beamB3.refractBeam('line', LensTwoEntryEquation, closestFlag=True, RI1 = 1, RI2 = LensRI)
beamB4 = beamClass.Beam(LensTwoPosB1, functions.lineAngle(refractedBeam), 'beamB4')

LensTwoPosB2 = beamB4.propagateBeam('circle', LensTwoExitEquation, closestFlag=True)
refractedBeam = beamB4.refractBeam('circle', LensTwoExitEquation, closestFlag=True,RI1 = LensRI, RI2 = 1)
beamB5 = beamClass.Beam(LensTwoPosB2, functions.lineAngle(refractedBeam), 'beamB5')

LensThreePosB1 = beamB5.propagateBeam('circle', LensThreeEntryEquation, closestFlag=True)
refractedBeam = beamB5.refractBeam('circle', LensThreeEntryEquation, closestFlag=True, RI1 = 1, RI2 = LensRI)
beamB6 = beamClass.Beam(LensThreePosB1, functions.lineAngle(refractedBeam), 'beamB6')

LensThreeExitEquation = [1, 0, -(sld + LensOneTwoDistance.x[0] + LensTwoThreeDistance)]

LensThreePosB2 = beamB6.propagateBeam('line', LensThreeExitEquation, closestFlag=True)
refractedBeam = beamB6.refractBeam('line', LensThreeExitEquation, closestFlag=True,RI1 = LensRI, RI2 = 1)
beamB7 = beamClass.Beam(LensThreePosB2, functions.lineAngle(refractedBeam), 'beamB7')


# =============== EQ-77 step 5 =============== #
""" propagate from L3 to mirror
"""

LensThreeMirrorDistance = 51.3       # variable

MirrorHypotenuse = 45
MirrorLength = np.sqrt(MirrorHypotenuse**2 / 2)
print(f'MirrorLength: {MirrorLength}')
MirrorTilt = np.radians(2)

def part4(variable):
    global beamA8, beamA9, beamA10, beamA11, beamA12
    global beamB8, beamB9, beamB10, beamB11, beamB12

    LensThreeMirrorDistance = variable[0]
    MirrorTilt = variable[1]

    MirrorPosition = [sld + LensOneTwoDistance.x[0] + LensTwoThreeDistance + LensThreeMirrorDistance, 0]
    MirrorEntryCorner1 = [ MirrorPosition[0] - np.sin(MirrorTilt)*MirrorLength/2, MirrorPosition[1] + np.cos(MirrorTilt)*MirrorLength/2]
    MirrorEntryCorner2 = [MirrorPosition[0] + np.sin(MirrorTilt)*MirrorLength/2, MirrorPosition[1] - np.cos(MirrorTilt)*MirrorLength/2]
    MirrorExitCorner1 = MirrorEntryCorner2
    MirrorExitCorner2 = [MirrorExitCorner1[0] + np.cos(MirrorTilt)*MirrorLength, MirrorExitCorner1[1] + np.sin(MirrorTilt)*MirrorLength]

    MirrorEntryLineEqn = functions.lineEquationPoints(MirrorEntryCorner1, MirrorEntryCorner2)
    MirrorExitLineEqn = functions.lineEquationPoints(MirrorExitCorner1, MirrorExitCorner2)
    MirrorHypotenusEqn = functions.lineEquationPoints(MirrorEntryCorner1, MirrorExitCorner2)

    """
    print(f'MirrorPosition: {MirrorPosition}')
    print(f'MirrorEntryCorner1: {MirrorEntryCorner1}')
    print(f'MirrorEntryCorner2: {MirrorEntryCorner2}')
    print(f'MirrorExitCorner1: {MirrorExitCorner1}')
    print(f'MirrorExitCorner2: {MirrorExitCorner2}')
    print(f'MirrorTilt: {MirrorTilt}')
    print(f'LensThreeMirrorDistance: {LensThreeMirrorDistance}')
    """
    # propagate to mirror
    MirrorPosA1 = beamA7.propagateBeam('line', MirrorEntryLineEqn, closestFlag=True)
    refractedBeam = beamA7.refractBeam('line', MirrorEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = 1)
    beamA8 = beamClass.Beam(MirrorPosA1, functions.lineAngle(refractedBeam), 'beamA8')
    MirrorPosA2 = beamA8.propagateBeam('line', MirrorHypotenusEqn, closestFlag=True)

    reflectedBeam = beamA8.reflectBeam('line', MirrorHypotenusEqn, closestFlag=True)
    beamA9 = beamClass.Beam(MirrorPosA2, functions.lineAngle(reflectedBeam), 'beamA9')
    MirrorPosA3 = beamA9.propagateBeam('line', MirrorExitLineEqn, closestFlag=True)
    refractedBeam = beamA9.refractBeam('line', MirrorExitLineEqn, closestFlag=True, RI1 = MirrorRI, RI2 = 1)
    beamA10 = beamClass.Beam(MirrorPosA3, functions.lineAngle(refractedBeam), 'beamA10')

    MirrorPosB1 = beamB7.propagateBeam('line', MirrorEntryLineEqn, closestFlag=True)
    refractedBeam = beamB7.refractBeam('line', MirrorEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = MirrorRI)
    beamB8 = beamClass.Beam(MirrorPosB1, functions.lineAngle(refractedBeam), 'beamB8')
    MirrorPosB2 = beamB8.propagateBeam('line', MirrorHypotenusEqn, closestFlag=True)

    reflectedBeam = beamB8.reflectBeam('line', MirrorHypotenusEqn, closestFlag=True)
    beamB9 = beamClass.Beam(MirrorPosB2, functions.lineAngle(reflectedBeam), 'beamB9')
    MirrorPosB3 = beamB9.propagateBeam('line', MirrorExitLineEqn, closestFlag=True)
    refractedBeam = beamB9.refractBeam('line', MirrorExitLineEqn, closestFlag=True, RI1 = MirrorRI, RI2 = 1)
    beamB10 = beamClass.Beam(MirrorPosB3, functions.lineAngle(refractedBeam), 'beamB10')
    # propagate to mirror

    # set up parameters: Mirror to sample
    #topMirrorWindowDistance = 234.7
    #MirrorMirrorDistance = 120                                                # can be varied in between optimization runs
    MirrorCentralDistance = 50                                               # referring to the Mirror in this simulation
    windowThickness = 2.997
    bottomWindowSampleDistance = 162.303

    # referring to the Mirror in this simulation
    MirrorWindowDistance = 188.10990
    centralLineOriginDistance = sld + LensOneTwoDistance.x[0] + LensTwoThreeDistance + LensThreeMirrorDistance + MirrorCentralDistance

    # referring to the Mirror in this simulation

    centralLineEqn = [1, 0, -(centralLineOriginDistance)]
    windowEntryLineEqn = [0, 1, MirrorWindowDistance]
    windowExitLineEqn = [0, 1, MirrorWindowDistance + windowThickness]
    targetHoriztonalLineEqn = [0, 1, (windowExitLineEqn[2]+bottomWindowSampleDistance)]
    targetPosition = functions.lineLineIntersection(centralLineEqn, targetHoriztonalLineEqn)
    # set up parameters: Mirror to sample
    
    # mirror to sample
    windowPosA1 = beamA10.propagateBeam('line', windowEntryLineEqn, closestFlag=True )
    refractedBeam = beamA10.refractBeam('line', windowEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = WindowsRI)
    beamA11 = beamClass.Beam(windowPosA1, functions.lineAngle(refractedBeam), 'beamA11')
    windowPosA2 = beamA11.propagateBeam('line', windowExitLineEqn, closestFlag=True)
    refractedBeam = beamA11.refractBeam('line', windowExitLineEqn, closestFlag=True, RI1 = WindowsRI, RI2 = 1)
    beamA12 = beamClass.Beam(windowPosA2, functions.lineAngle(refractedBeam), 'beamA12')
    samplePositionA = beamA12.propagateBeam('line', targetHoriztonalLineEqn, closestFlag=True)
    
    windowPosB1 = beamB10.propagateBeam('line', windowEntryLineEqn, closestFlag=True )
    refractedBeam = beamB10.refractBeam('line', windowEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = WindowsRI)
    beamB11 = beamClass.Beam(windowPosB1, functions.lineAngle(refractedBeam), 'beamB11')
    windowPosB2 = beamB11.propagateBeam('line', windowExitLineEqn, closestFlag=True)
    refractedBeam = beamB11.refractBeam('line', windowExitLineEqn, closestFlag=True, RI1 = WindowsRI, RI2 = 1)
    beamB12 = beamClass.Beam(windowPosB2, functions.lineAngle(refractedBeam), 'beamB12')
    samplePositionB = beamB12.propagateBeam('line', targetHoriztonalLineEqn, closestFlag=True)
    # mirror to sample

    targetSampleDistanceA = samplePositionA[0] - targetPosition[0]
    targetSampleDistanceB = samplePositionB[0] - targetPosition[0]
    beamABDistance = samplePositionA[0] - samplePositionB[0]

    return (np.abs(beamABDistance) + np.abs(targetSampleDistanceA))

LensThreeMirrorDistance = minimize(part4, (200,0.035), tol=1e-10, method = 'Nelder-Mead')

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


beamB1.printBeam()
beamB2.printBeam()
beamB3.printBeam()
beamB4.printBeam()
beamB5.printBeam()
beamB6.printBeam()
beamB7.printBeam()
beamB8.printBeam()
beamB9.printBeam()
beamB10.printBeam()
beamB11.printBeam()
beamB12.printBeam()

print(LensThreeMirrorDistance)