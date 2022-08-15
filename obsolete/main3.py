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
# set up parameters: E2 lens
LensTwoFullHeight = 25.4
LensTwoRadius = 21.7
LensTwoEdgeThickness = 2
LensTwoCentreThickness = functions.calculatePCLensCT(LensTwoFullHeight, LensTwoRadius, LensTwoEdgeThickness)

def step2part2(LensOneTwoDistance):

    global beam3, beam4, beam5

    LensOneTwoDistance = LensOneTwoDistance[0]
    LensTwoEntryEquation = [1,0,LensOneTwoDistance+SourcesLensDistance.x[0]]
    LensTwoPos1 = beam3.propagateBeam('line', LensTwoEntryEquation, closestFlag=True)
    refractedBeam = beam3.refractBeam('line', LensTwoEntryEquation, closestFlag=True, RI1 = 1, RI2 = LensRI)
    beam4 = beamClass.Beam(LensTwoPos1, functions.lineAngle(refractedBeam), 'beam4')

    LensTwoExitEquation = [-(SourcesLensDistance.x[0] + LensOneTwoDistance - (LensTwoRadius - LensTwoCentreThickness)), 0, LensTwoRadius]

    LensTwoPos2 = beam4.propagateBeam('circle', LensTwoExitEquation, closestFlag=True)
    refractedBeam = beam4.refractBeam('circle', LensTwoExitEquation, closestFlag=True,RI1 = LensRI, RI2 = 1)
    beam5 = beamClass.Beam(LensTwoPos2, functions.lineAngle(refractedBeam), 'beam5')

    return np.abs(beam5.gradient)

LensOneTwoDistance = minimize(step2part2, 2, tol=1e-10, method = 'Nelder-Mead')

beam3.printBeam()
beam4.printBeam()
print(LensOneTwoDistance)


# =============== EQ-step2-part3 =============== #
# set up parameters: E3 lens
LensThreeFullHeight = 25.4
LensThreeRadius = 216.9
LensThreeEdgeThickness = 2
LensThreeCentreThickness = functions.calculatePCLensCT(LensThreeFullHeight, LensThreeRadius, LensThreeEdgeThickness)

LensTwoThreeDistance = 50
LensThreeEntryEquation = [-(SourcesLensDistance.x[0] + LensOneTwoDistance.x[0] + LensTwoThreeDistance + LensThreeRadius - LensThreeCentreThickness), 0, LensThreeRadius]

LensThreePos1 = beam5.propagateBeam('circle', LensThreeEntryEquation, closestFlag=True)
refractedBeam = beam5.refractBeam('circle', LensThreeEntryEquation, closestFlag=True, RI1 = 1, RI2 = LensRI)
beam6 = beamClass.Beam(LensThreePos1, functions.lineAngle(refractedBeam), 'beam6')

LensThreeExitEquation = [1, 0, (SourcesLensDistance.x[0] + LensOneTwoDistance.x[0] + LensTwoThreeDistance)]

LensThreePos2 = beam6.propagateBeam('line', LensThreeExitEquation, closestFlag=True)
refractedBeam = beam6.refractBeam('line', LensThreeExitEquation, closestFlag=True,RI1 = LensRI, RI2 = 1)
beam7 = beamClass.Beam(LensThreePos2, functions.lineAngle(refractedBeam), 'beam7')

beam5.printBeam()
beam6.printBeam()

def step2part3(angle):
    
    global beam7, beam8, beam9, beam10, beam11, beam12

    angle = angle[0]

    # set up parameters: prism related
    PrismLength = 25.4
    PrismTilt = np.radians(angle)                           # variable

    LensThreePrismDistance = 50

    PrismPosition = [-(SourcesLensDistance.x[0] + LensOneTwoDistance.x[0] + LensTwoThreeDistance + LensThreePrismDistance), 0]
    PrismEntryCorner1 = [PrismPosition[0] - np.sin(PrismTilt)*PrismLength/2, -(PrismPosition[1] + np.cos(PrismTilt)*PrismLength/2)]
    PrismEntryCorner2 = [PrismPosition[0] + np.sin(PrismTilt)*PrismLength/2, -(PrismPosition[1] - np.cos(PrismTilt)*PrismLength/2)]
    PrismExitCorner1 = PrismEntryCorner1
    PrismExitCorner2 = [-(np.abs(PrismExitCorner1[0]) + np.cos(PrismTilt)*PrismLength), -(np.abs(PrismExitCorner1[1]) - np.sin(PrismTilt)*PrismLength) ]

    PrismEntryLineEqn = functions.lineEquationPoints(PrismEntryCorner1, PrismEntryCorner2)
    PrismExitLineEqn = functions.lineEquationPoints(PrismExitCorner1, PrismExitCorner2)
    PrismHypotenusEqn = functions.lineEquationPoints(PrismEntryCorner2, PrismExitCorner2)

    PrismPos1 = beam7.propagateBeam('line', PrismEntryLineEqn, closestFlag=True)
    refractedBeam = beam7.refractBeam('line', PrismEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = PrismRI)
    beam8 = beamClass.Beam(PrismPos1, functions.lineAngle(refractedBeam), 'beam8')

    PrismPos2 = beam8.propagateBeam('line', PrismHypotenusEqn, closestFlag=True)
    reflectedBeam = beam8.reflectBeam('line', PrismHypotenusEqn, closestFlag=True)
    beam9 = beamClass.Beam(PrismPos2, functions.lineAngle(reflectedBeam), 'beam9')

    PrismPos3 = beam9.propagateBeam('line', PrismExitLineEqn, closestFlag=True)
    refractedBeam = beam9.refractBeam('line', PrismExitLineEqn, closestFlag=True, RI1 = PrismRI, RI2 = 1)
    beam10 = beamClass.Beam(PrismPos3, functions.lineAngle(refractedBeam), 'beam10')

    # set up parameters: prism to sample
    topPrismWindowDistance = 234.7
    prismPrismDistance = 65                                                  # can be varied in between optimization runs
    prismCentralDistance = 45                                                # referring to the prism in this simulation
    windowThickness = 2.997
    bottomWindowSampleDistance = 231.703

    # referring to the prism in this simulation
    prismWindowDistance = topPrismWindowDistance - prismPrismDistance
    centralLineOriginDistance = np.abs(PrismPosition[0]) + 45
    # referring to the prism in this simulation

    centralLineEqn = [1, 0, centralLineOriginDistance]
    windowEntryLineEqn = [0, 1, prismWindowDistance]
    windowExitLineEqn = [0, 1, prismWindowDistance + windowThickness]
    targetHoriztonalLineEqn = [0, 1, (windowExitLineEqn[2] + bottomWindowSampleDistance)]

    targetPosition = functions.lineLineIntersection(centralLineEqn, targetHoriztonalLineEqn)

    windowPos1 = beam10.propagateBeam('line', windowEntryLineEqn, closestFlag=True )
    refractedBeam = beam10.refractBeam('line', windowEntryLineEqn, closestFlag=True, RI1 = 1, RI2 = WindowsRI)
    beam11 = beamClass.Beam(windowPos1, functions.lineAngle(refractedBeam), 'beam11')
    windowPos2 = beam11.propagateBeam('line', windowExitLineEqn, closestFlag=True)
    refractedBeam = beam11.refractBeam('line', windowExitLineEqn, closestFlag=True, RI1 = WindowsRI, RI2 = 1)
    beam12 = beamClass.Beam(windowPos2, functions.lineAngle(refractedBeam), 'beam12')
    samplePosition = beam12.propagateBeam('line', targetHoriztonalLineEqn, closestFlag=True)

    targetSampleDistance = samplePosition[0] - targetPosition[0]

    return np.abs(targetSampleDistance)

angle = minimize(step2part3, 1, tol=1e-10, method = 'Nelder-Mead')

beam7.printBeam()
beam8.printBeam()
beam9.printBeam()
beam10.printBeam()
beam11.printBeam()
beam12.printBeam()

print(angle)