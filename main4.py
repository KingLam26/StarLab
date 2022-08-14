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
    global beamA1, beamA2, beamA3, sld

    sld = SourcesLensDistance[0]

    beamA1 = beamClass.Beam([0, 0], np.radians(14), 'beamA1')
    LensOnePosA1 = beamA1.propagateBeam('line', [1,0,-sld], closestFlag=True)
    beamA1.refractBeam('line', [1,0,-sld], closestFlag=True, RI1 = 1, RI2 = LensRI)

    beamA2 = beamClass.Beam(LensOnePosA1, beamA1.ExitAngle, 'beamA2')
    LensOnePosA2 = beamA2.propagateBeam('circle', [sld - (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True)
    refracted = beamA2.refractBeam('circle', [sld - (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True, RI1 = LensRI, RI2 = 1)

    beamA3 = beamClass.Beam(LensOnePosA2, functions.lineAngle(refracted), 'beamA3')
    return np.abs(beamA3.gradient-0)

SourcesLensDistance = minimize(step2part1, 2, tol=1e-10, method = 'Nelder-Mead')

# print information about the beams
print(f'SourcesLensDistance: {SourcesLensDistance.x[0]}')


# =============== Horiba-step2-part2 =============== #
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


# =============== L12542-step2-part3.1 =============== #
# it is at this point, that we need to build the other beam, up to beam B4...

beamB1 = beamClass.Beam([0, 0], np.radians(-14), 'beamB1')
LensOnePosB1 = beamB1.propagateBeam('line', [1,0,-sld], closestFlag=True)
beamB1.refractBeam('line', [1,0,-sld], closestFlag=True, RI1 = 1, RI2 = LensRI)

beamB2 = beamClass.Beam(LensOnePosB1, -beamA1.ExitAngle, 'beamB2')
LensOnePosB2 = beamB2.propagateBeam('circle', [sld - (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True)
refracted = beamA2.refractBeam('circle', [sld - (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True, RI1 = LensRI, RI2 = 1)

beamB3 = beamClass.Beam(LensOnePosB2, functions.lineAngle(refracted), 'beamB3')
LensTwoPosB1 = beamB3.propagateBeam('circle', LensTwoEntryEquation, closestFlag=True)
refractedBeam = beamB3.refractBeam('circle', LensTwoEntryEquation, closestFlag=True, RI1 = 1, RI2 = LensRI)
beamB4 = beamClass.Beam(LensTwoPosB1, functions.lineAngle(refractedBeam), 'beamB4')

LensTwoPosB2 = beamB4.propagateBeam('line', [1,0,-(SourcesLensDistance.x[0] + LensOneTwoDistance)], closestFlag=True)
refractedBeam = beamB4.refractBeam('line', [1,0,-(SourcesLensDistance.x[0] + LensOneTwoDistance)], closestFlag=True,RI1 = LensRI, RI2 = 1)
beamB5 = beamClass.Beam(LensTwoPosB2, functions.lineAngle(refractedBeam), 'beamB5')


# =============== L12542-step2-part3.2 =============== #
# set up parameters: prism related
PrismLength = 25.4
PrismTilt = 0.04698284
LensTwoPrismDistance = 50                   # variable

"""
PrismPosition = [SourcesLensDistance.x[0] + LensOneTwoDistance + LensTwoPrismDistance, 0]
PrismEntryCorner1 = [PrismPosition[0] - np.sin(PrismTilt)*PrismLength/2, PrismPosition[1] + np.cos(PrismTilt)*PrismLength/2]
PrismEntryCorner2 = [PrismPosition[0] + np.sin(PrismTilt)*PrismLength/2, PrismPosition[1] - np.cos(PrismTilt)*PrismLength/2]
PrismExitCorner1 = PrismEntryCorner2
PrismExitCorner2 = [PrismExitCorner1[0] + np.cos(PrismTilt)*PrismLength, -(np.abs(PrismExitCorner1[1]) - np.sin(PrismTilt)*PrismLength) ]

PrismEntryLineEqn = functions.lineEquationPoints(PrismEntryCorner1, PrismEntryCorner2)
PrismExitLineEqn = functions.lineEquationPoints(PrismExitCorner1, PrismExitCorner2)
PrismHypotenusEqn = functions.lineEquationPoints(PrismEntryCorner1, PrismExitCorner2)

PrismPosB1 = beamB5.propagateBeam('line', PrismEntryLineEqn, closestFlag=True)
refractedBeam = beamB5.refractBeam('line', PrismEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = PrismRI)
beamB6 = beamClass.Beam(PrismPosB1, functions.lineAngle(refractedBeam), 'beamB6')

PrismPosB2 = beamB6.propagateBeam('line', PrismHypotenusEqn, closestFlag=True)
reflectedBeam = beamB6.reflectBeam('line', PrismHypotenusEqn, closestFlag=True)
beamB7 = beamClass.Beam(PrismPosB2, functions.lineAngle(reflectedBeam), 'beamB7')

PrismPosB3 = beamB7.propagateBeam('line', PrismExitLineEqn, closestFlag=True)
refractedBeam = beamB7.refractBeam('line', PrismExitLineEqn, closestFlag=True, RI1 = PrismRI, RI2 = 1)
beamB8 = beamClass.Beam(PrismPosB3, functions.lineAngle(refractedBeam), 'beamB8')

# set up parameters: prism to sample
topPrismWindowDistance = 234.7
prismPrismDistance = 120                                                # can be varied in between optimization runs
prismCentralDistance = 45                                               # referring to the prism in this simulation
windowThickness = 2.997
bottomWindowSampleDistance = 231.703

# referring to the prism in this simulation
prismWindowDistance = topPrismWindowDistance - prismPrismDistance
centralLineOriginDistance = SourcesLensDistance.x[0] + LensOneTwoDistance + LensTwoPrismDistance + prismCentralDistance
# referring to the prism in this simulation

centralLineEqn = [1, 0, -(centralLineOriginDistance)]
windowEntryLineEqn = [0, 1, prismWindowDistance]
windowExitLineEqn = [0, 1, prismWindowDistance + windowThickness]
targetHoriztonalLineEqn = [0, 1, (windowExitLineEqn[2]+bottomWindowSampleDistance)]

targetPosition = functions.lineLineIntersection(centralLineEqn, targetHoriztonalLineEqn)

windowPosB1 = beamB8.propagateBeam('line', windowEntryLineEqn, closestFlag=True )
refractedBeam = beamB8.refractBeam('line', windowEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = WindowsRI)
beamB9 = beamClass.Beam(windowPosB1, functions.lineAngle(refractedBeam), 'beamB9')
windowPosB2 = beamB9.propagateBeam('line', windowExitLineEqn, closestFlag=True)
refractedBeam = beamB9.refractBeam('line', windowExitLineEqn, closestFlag=True, RI1 = WindowsRI, RI2 = 1)
beamB10 = beamClass.Beam(windowPosB2, functions.lineAngle(refractedBeam), 'beamA10')
samplePosB = beamB10.propagateBeam('line', targetHoriztonalLineEqn, closestFlag=True)

beamB4.printBeam()
beamB5.printBeam()
beamB6.printBeam()
beamB7.printBeam()
beamB8.printBeam()
beamB9.printBeam()
beamB10.printBeam()
"""


# define function
def step2part3(var):

    global beamA5, beamA6, beamA7, beamA8, beamA9, beamA10
    global beamB5, beamB6, beamB7, beamB8, beamB9, beamB10

    global targetSampleDistanceA, targetSampleDistanceB, beamABDistance

    LensTwoPrismDistance = var[0]
    PrismTilt = var[1]

    # calculating prism parameters
    PrismPosition = [SourcesLensDistance.x[0] + LensOneTwoDistance + LensTwoPrismDistance, 0]
    PrismEntryCorner1 = [PrismPosition[0] - np.sin(PrismTilt)*PrismLength/2, PrismPosition[1] + np.cos(PrismTilt)*PrismLength/2]
    PrismEntryCorner2 = [PrismPosition[0] + np.sin(PrismTilt)*PrismLength/2, PrismPosition[1] - np.cos(PrismTilt)*PrismLength/2]
    PrismExitCorner1 = PrismEntryCorner2
    PrismExitCorner2 = [PrismExitCorner1[0] + np.cos(PrismTilt)*PrismLength, -(np.abs(PrismExitCorner1[1]) - np.sin(PrismTilt)*PrismLength) ]

    PrismEntryLineEqn = functions.lineEquationPoints(PrismEntryCorner1, PrismEntryCorner2)
    PrismExitLineEqn = functions.lineEquationPoints(PrismExitCorner1, PrismExitCorner2)
    PrismHypotenusEqn = functions.lineEquationPoints(PrismEntryCorner1, PrismExitCorner2)
    # calculating prism parameters

    # handling beamA at the prism
    PrismPosA1 = beamA5.propagateBeam('line', PrismEntryLineEqn, closestFlag=True)

    refractedBeam = beamA5.refractBeam('line', PrismEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = PrismRI)
    beamA6 = beamClass.Beam(PrismPosA1, functions.lineAngle(refractedBeam), 'beamA6')
    PrismPosA2 = beamA6.propagateBeam('line', PrismHypotenusEqn, closestFlag=True)

    reflectedBeam = beamA6.reflectBeam('line', PrismHypotenusEqn, closestFlag=True)
    beamA7 = beamClass.Beam(PrismPosA2, functions.lineAngle(reflectedBeam), 'beamA7')
    PrismPosA3 = beamA7.propagateBeam('line', PrismExitLineEqn, closestFlag=True)

    refractedBeam = beamA7.refractBeam('line', PrismExitLineEqn, closestFlag=True, RI1 = PrismRI, RI2 = 1)
    beamA8 = beamClass.Beam(PrismPosA3, functions.lineAngle(refractedBeam), 'beamA8')
    # handling beamA at the prism

    # handling beamB at the prism
    PrismPosB1 = beamB5.propagateBeam('line', PrismEntryLineEqn, closestFlag=True)
    refractedBeam = beamB5.refractBeam('line', PrismEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = PrismRI)
    beamB6 = beamClass.Beam(PrismPosB1, functions.lineAngle(refractedBeam), 'beamB6')

    PrismPosB2 = beamB6.propagateBeam('line', PrismHypotenusEqn, closestFlag=True)
    reflectedBeam = beamB6.reflectBeam('line', PrismHypotenusEqn, closestFlag=True)
    beamB7 = beamClass.Beam(PrismPosB2, functions.lineAngle(reflectedBeam), 'beamB7')

    PrismPosB3 = beamB7.propagateBeam('line', PrismExitLineEqn, closestFlag=True)
    refractedBeam = beamB7.refractBeam('line', PrismExitLineEqn, closestFlag=True, RI1 = PrismRI, RI2 = 1)
    beamB8 = beamClass.Beam(PrismPosB3, functions.lineAngle(refractedBeam), 'beamB8')
    # handling beamB at the prism

    # set up parameters: prism to sample
    topPrismWindowDistance = 234.7
    prismPrismDistance = 120                                                # can be varied in between optimization runs
    prismCentralDistance = 45                                               # referring to the prism in this simulation
    windowThickness = 2.997
    bottomWindowSampleDistance = 231.703

    # referring to the prism in this simulation
    prismWindowDistance = topPrismWindowDistance - prismPrismDistance
    centralLineOriginDistance = SourcesLensDistance.x[0] + LensOneTwoDistance + LensTwoPrismDistance + prismCentralDistance
    # referring to the prism in this simulation

    centralLineEqn = [1, 0, -(centralLineOriginDistance)]
    windowEntryLineEqn = [0, 1, prismWindowDistance]
    windowExitLineEqn = [0, 1, prismWindowDistance + windowThickness]
    targetHoriztonalLineEqn = [0, 1, (windowExitLineEqn[2]+bottomWindowSampleDistance)]
    targetPosition = functions.lineLineIntersection(centralLineEqn, targetHoriztonalLineEqn)
    # set up parameters: prism to sample


    # handling beamA, prism to sample
    windowPosA1 = beamA8.propagateBeam('line', windowEntryLineEqn, closestFlag=True )
    refractedBeam = beamA8.refractBeam('line', windowEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = WindowsRI)
    beamA9 = beamClass.Beam(windowPosA1, functions.lineAngle(refractedBeam), 'beamA9')
    windowPosA2 = beamA9.propagateBeam('line', windowExitLineEqn, closestFlag=True)
    refractedBeam = beamA9.refractBeam('line', windowExitLineEqn, closestFlag=True, RI1 = WindowsRI, RI2 = 1)
    beamA10 = beamClass.Beam(windowPosA2, functions.lineAngle(refractedBeam), 'beamA10')
    samplePosA = beamA10.propagateBeam('line', targetHoriztonalLineEqn, closestFlag=True)
    # handling beamA, prism to sample

    # handling beamB, prism to sample
    windowPosB1 = beamB8.propagateBeam('line', windowEntryLineEqn, closestFlag=True )
    refractedBeam = beamB8.refractBeam('line', windowEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = WindowsRI)
    beamB9 = beamClass.Beam(windowPosB1, functions.lineAngle(refractedBeam), 'beamB9')
    windowPosB2 = beamB9.propagateBeam('line', windowExitLineEqn, closestFlag=True)
    refractedBeam = beamB9.refractBeam('line', windowExitLineEqn, closestFlag=True, RI1 = WindowsRI, RI2 = 1)
    beamB10 = beamClass.Beam(windowPosB2, functions.lineAngle(refractedBeam), 'beamB10')
    samplePosB = beamB10.propagateBeam('line', targetHoriztonalLineEqn, closestFlag=True)
    # handling beamB, prism to sample

    # calculating distance between beamA and beamB on the horizontal sample line
    targetSampleDistanceA = samplePosA[0] - targetPosition[0]
    targetSampleDistanceB = samplePosB[0] - targetPosition[0]
    beamABDistance = samplePosA[0] - samplePosB[0]
    # calculating distance between beamA and beamB on the horizontal sample line

    return (np.abs(beamABDistance) + np.abs(targetSampleDistanceA))

LensTwoPrismDistance = minimize(step2part3, (40,0.01), tol=1e-10, method = 'Nelder-Mead')

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

print(f'targetSampleDistanceA: {targetSampleDistanceA}')
print(f'beamABDistance: {beamABDistance}')

print(f'PrismTilt: {LensTwoPrismDistance}')