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
# set up parameters: E1 lens
sysTwoOrigin = []
LensOneFullHeight = 12.7
LensOneRadius = 8.7
LensOneEdgeThickness = 1.5
LensOneCentreThickness = functions.calculatePCLensCT(LensOneFullHeight, LensOneRadius, LensOneEdgeThickness)
clearApertureRadius = 0.9*LensOneFullHeight/2

# define function
def step2part1(SourcesLensDistance):
    global beamA1, beamA2, beamA3

    global LensOnePosA1, LensOnePosA2, sld

    sld = SourcesLensDistance[0]

    beamA1 = beamClass.Beam([0, 0], np.radians(150), 'beamA1')
    LensOnePosA1 = beamA1.propagateBeam('line', [1,0,sld], closestFlag=True)
    beamA1.refractBeam('line', [1,0,sld], closestFlag=True, RI1 = 1, RI2 = LensRI)

    beamA2 = beamClass.Beam(LensOnePosA1, np.pi - beamA1.ExitAngle, 'beamA2')
    LensOnePosA2 = beamA2.propagateBeam('circle', [-sld + (LensOneRadius - LensOneCentreThickness), 0, LensOneRadius], closestFlag=True)
    refracted = beamA2.refractBeam('circle', [-sld + (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True, RI1 = LensRI, RI2 = 1)

    beamA3 = beamClass.Beam(LensOnePosA2, functions.lineAngle(refracted), 'beamA3')

    return np.abs(LensOnePosA2[1] - clearApertureRadius)

SourcesLensDistance = minimize(step2part1, 2, tol=1e-10, method = 'Nelder-Mead')

# print information about the beams
print(f'SourcesLensDistance: {SourcesLensDistance.x[0]}')


# =============== EQ-step2-part2 =============== #
# set up parameters: E2 lens
LensTwoFullHeight = 25.4
LensTwoRadius = 21.7
LensTwoEdgeThickness = 2
LensTwoCentreThickness = functions.calculatePCLensCT(LensTwoFullHeight, LensTwoRadius, LensTwoEdgeThickness)

def step2part2(LensOneTwoDistance):

    global beamA3, beamA4, beamA5

    global LensTwoPosA1, LensTwoPosA2

    global LensTwoEntryEquation, LensTwoExitEquation

    LensOneTwoDistance = LensOneTwoDistance[0]
    LensTwoEntryEquation = [1,0,LensOneTwoDistance+SourcesLensDistance.x[0]]
    LensTwoPosA1 = beamA3.propagateBeam('line', LensTwoEntryEquation, closestFlag=True)
    refractedBeam = beamA3.refractBeam('line', LensTwoEntryEquation, closestFlag=True, RI1 = 1, RI2 = LensRI)
    beamA4 = beamClass.Beam(LensTwoPosA1, functions.lineAngle(refractedBeam), 'beamA4')

    LensTwoExitEquation = [-(SourcesLensDistance.x[0] + LensOneTwoDistance - (LensTwoRadius - LensTwoCentreThickness)), 0, LensTwoRadius]

    LensTwoPosA2 = beamA4.propagateBeam('circle', LensTwoExitEquation, closestFlag=True)
    refractedBeam = beamA4.refractBeam('circle', LensTwoExitEquation, closestFlag=True,RI1 = LensRI, RI2 = 1)
    beamA5 = beamClass.Beam(LensTwoPosA2, functions.lineAngle(refractedBeam), 'beamA5')

    return np.abs(beamA5.gradient)

LensOneTwoDistance = minimize(step2part2, 2, tol=1e-10, method = 'Nelder-Mead')

print(LensOneTwoDistance)


# =============== EQ-step2-part3 =============== #
# set up parameters: E3 lens
LensThreeFullHeight = 25.4
LensThreeRadius = 216.9
LensThreeEdgeThickness = 2
LensThreeCentreThickness = functions.calculatePCLensCT(LensThreeFullHeight, LensThreeRadius, LensThreeEdgeThickness)

LensTwoThreeDistance = 50
LensThreeEntryEquation = [-(SourcesLensDistance.x[0] + LensOneTwoDistance.x[0] + LensTwoThreeDistance + LensThreeRadius - LensThreeCentreThickness), 0, LensThreeRadius]

LensThreePosA1 = beamA5.propagateBeam('circle', LensThreeEntryEquation, closestFlag=True)
refractedBeam = beamA5.refractBeam('circle', LensThreeEntryEquation, closestFlag=True, RI1 = 1, RI2 = LensRI)
beamA6 = beamClass.Beam(LensThreePosA1, functions.lineAngle(refractedBeam), 'beamA6')

LensThreeExitEquation = [1, 0, (SourcesLensDistance.x[0] + LensOneTwoDistance.x[0] + LensTwoThreeDistance)]

LensThreePosA2 = beamA6.propagateBeam('line', LensThreeExitEquation, closestFlag=True)
refractedBeam = beamA6.refractBeam('line', LensThreeExitEquation, closestFlag=True,RI1 = LensRI, RI2 = 1)
beamA7 = beamClass.Beam(LensThreePosA2, functions.lineAngle(refractedBeam), 'beamA7')


# =============== EQ-step2-part4.1 =============== #
# here we construct the other beam, up to beam B7
beamB1 = beamClass.Beam([0, 0], np.radians(-150), 'beamB1')
LensOnePosB1 = beamB1.propagateBeam('line', [1,0,sld], closestFlag=True)
beamB1.refractBeam('line', [1,0,sld], closestFlag=True, RI1 = 1, RI2 = LensRI)

beamB2 = beamClass.Beam(LensOnePosB1, -(np.pi - beamA1.ExitAngle), 'beamB2')
LensOnePosB2 = beamB2.propagateBeam('circle', [-sld + (LensOneRadius - LensOneCentreThickness), 0, LensOneRadius], closestFlag=True)
refracted = beamB2.refractBeam('circle', [-sld + (LensOneRadius - LensOneCentreThickness),0,LensOneRadius], closestFlag=True, RI1 = LensRI, RI2 = 1)
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

LensThreePosB2 = beamB6.propagateBeam('line', LensThreeExitEquation, closestFlag=True)
refractedBeam = beamB6.refractBeam('line', LensThreeExitEquation, closestFlag=True,RI1 = LensRI, RI2 = 1)
beamB7 = beamClass.Beam(LensThreePosB2, functions.lineAngle(refractedBeam), 'beamB7')

LensThreePrismDistance = 127                 # variable 2
PrismTilt = np.radians(2.55)

# set up parameters: prism related
PrismLength = 25.4

PrismPosition = [-(SourcesLensDistance.x[0] + LensOneTwoDistance.x[0] + LensTwoThreeDistance + LensThreePrismDistance), 0]
PrismEntryCorner1 = [PrismPosition[0] - np.sin(PrismTilt)*PrismLength/2, -(PrismPosition[1] + np.cos(PrismTilt)*PrismLength/2)]
PrismEntryCorner2 = [PrismPosition[0] + np.sin(PrismTilt)*PrismLength/2, -(PrismPosition[1] - np.cos(PrismTilt)*PrismLength/2)]
PrismExitCorner1 = PrismEntryCorner1
PrismExitCorner2 = [-(np.abs(PrismExitCorner1[0]) + np.cos(PrismTilt)*PrismLength), -(np.abs(PrismExitCorner1[1]) - np.sin(PrismTilt)*PrismLength) ]

PrismEntryLineEqn = functions.lineEquationPoints(PrismEntryCorner1, PrismEntryCorner2)
PrismExitLineEqn = functions.lineEquationPoints(PrismExitCorner1, PrismExitCorner2)
PrismHypotenusEqn = functions.lineEquationPoints(PrismEntryCorner2, PrismExitCorner2)
# set up parameters: prism related

# propogate beam A: E3 to prism
prismPosA1 = beamA7.propagateBeam('line', PrismEntryLineEqn, closestFlag=True)
refractedBeam = beamA7.refractBeam('line', PrismEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = PrismRI)
beamA8 = beamClass.Beam(prismPosA1, functions.lineAngle(refractedBeam), 'beamA8')

prismPosA2 = beamA8.propagateBeam('line', PrismHypotenusEqn, closestFlag=True)
reflectedBeam = beamA8.reflectBeam('line', PrismHypotenusEqn, closestFlag=True)
beamA9 = beamClass.Beam(prismPosA2, functions.lineAngle(reflectedBeam), 'beamA9')

prismPosA3 = beamA9.propagateBeam('line', PrismExitLineEqn, closestFlag=True)
refractedBeam = beamA9.refractBeam('line', PrismExitLineEqn, closestFlag=True, RI1 = PrismRI, RI2 = 1)
beamA10 = beamClass.Beam(prismPosA3, functions.lineAngle(refractedBeam), 'beamA10')
# propogate beam A: E3 to prism

# propogate beam B: E3 to prism
prismPosB1 = beamB7.propagateBeam('line', PrismEntryLineEqn, closestFlag=True)
refractedBeam = beamB7.refractBeam('line', PrismEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = PrismRI)
beamB8 = beamClass.Beam(prismPosB1, functions.lineAngle(refractedBeam), 'beamB8')

prismPosB2 = beamB8.propagateBeam('line', PrismHypotenusEqn, closestFlag=True)
reflectedBeam = beamB8.reflectBeam('line', PrismHypotenusEqn, closestFlag=True)
beamB9 = beamClass.Beam(prismPosB2, functions.lineAngle(reflectedBeam), 'beamB9')

prismPosB3 = beamB9.propagateBeam('line', PrismExitLineEqn, closestFlag=True)
refractedBeam = beamB9.refractBeam('line', PrismExitLineEqn, closestFlag=True, RI1 = PrismRI, RI2 = 1)
beamB10 = beamClass.Beam(prismPosB3, functions.lineAngle(refractedBeam), 'beamB10')
# propogate beam B: E3 to prism

# set up parameters: prism to sample and target
topPrismWindowDistance = 264.7
prismPrismDistance = 65                                                  # can be varied in between optimization runs, between lower surface of top prism
prismCentralDistance = 45                                                # referring to the prism in this simulation
windowThickness = 2.997
bottomWindowSampleDistance = 162.303

# referring to the prism in this simulation
# prismWindowDistance = topPrismWindowDistance - prismPrismDistance
prismWindowDistance = 199.7
prismCentralDistance = 45
centralLineOriginDistance = np.abs(PrismPosition[0]) + prismCentralDistance
# referring to the prism in this simulation

centralLineEqn = [1, 0, centralLineOriginDistance]
windowEntryLineEqn = [0, 1, prismWindowDistance]
windowExitLineEqn = [0, 1, prismWindowDistance + windowThickness]
targetHoriztonalLineEqn = [0, 1, (windowExitLineEqn[2] + bottomWindowSampleDistance)]

targetPosition = functions.lineLineIntersection(centralLineEqn, targetHoriztonalLineEqn)
# set up parameters: prism to sample and target

# propogate beam A: prism to sample
windowPosA1 = beamA10.propagateBeam('line', windowEntryLineEqn, closestFlag=True )
refractedBeam = beamA10.refractBeam('line', windowEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = WindowsRI)
beamA11 = beamClass.Beam(windowPosA1, functions.lineAngle(refractedBeam), 'beamA11')
windowPosA2 = beamA11.propagateBeam('line', windowExitLineEqn, closestFlag=True)
refractedBeam = beamA11.refractBeam('line', windowExitLineEqn, closestFlag=True, RI1 = WindowsRI, RI2 = 1)
beamA12 = beamClass.Beam(windowPosA2, functions.lineAngle(refractedBeam), 'beamA12')
samplePosA = beamA12.propagateBeam('line', targetHoriztonalLineEqn, closestFlag=True)
# propogate beam A: prism to sample

# propogate beam B: prism to sample
windowPosB1 = beamB10.propagateBeam('line', windowEntryLineEqn, closestFlag=True )
refractedBeam = beamB10.refractBeam('line', windowEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = WindowsRI)
beamB11 = beamClass.Beam(windowPosB1, functions.lineAngle(refractedBeam), 'beamB11')
windowPosB2 = beamB11.propagateBeam('line', windowExitLineEqn, closestFlag=True)
refractedBeam = beamB11.refractBeam('line', windowExitLineEqn, closestFlag=True, RI1 = WindowsRI, RI2 = 1)
beamB12 = beamClass.Beam(windowPosB2, functions.lineAngle(refractedBeam), 'beamB12')
samplePosB = beamB12.propagateBeam('line', targetHoriztonalLineEqn, closestFlag=True)
# propogate beam B: prism to sample

# handling the target and sample positions
targetSampleDistanceA = samplePosA[0] - targetPosition[0]
targetSampleDistanceB = samplePosB[0] - targetPosition[0]
beamABDistance = samplePosA[0] - samplePosB[0]
# handling the target and sample positions



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


print(samplePosA)
print(samplePosB)

print(f'beamABDistance: {beamABDistance}')
print(f'targetSampleDistanceA: {targetSampleDistanceA}')