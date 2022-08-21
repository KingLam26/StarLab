# =============== import modules =============== #
from math import radians
from tkinter import E
import numpy as np
from functions import *


# =============== Beam class =============== #

class Beam:
    def __init__(self, initialPos, theta, name):
        """ A beam is an infinite straight line defined by a pair of x-y coordinates and angle theta

        Args:
            name (str):                         name of the beam
            initialPos (list of floats):        [x,y] planar coordinates of a point on the line, usually the point where the beam begins
            theta (float):                      relative to the +ve x-axis, positive indicates counter-clockwise, between 0 and np.pi, radians
        """
        self._name = name
        self._initialPos = initialPos
        self.theta = theta
        self._gradient = self.beamGradient()
        self._intercept = self.beamintercept()
        self._beamEquation = lineEquationMC(self.gradient, self.intercept)     # [A,B,C]
        self.incidentAngle = np.nan
        self.ExitAngle = np.nan

    def beamGradient(self):
        while self.theta > np.pi:
            self.theta -= np.pi

        if self.theta == np.pi / 2:
            return np.inf
        else:
            return np.tan(self.theta)
    
    def beamintercept(self):
        if self.gradient == np.inf:
            return self.initialPos[0]
        else:
            return self.initialPos[1] - self.gradient*self.initialPos[0]

    def propagateBeam(self, surfaceType, surfaceEquation, closestFlag):
        
        """ propagates beam until it intersects a surface, either a line or circle
        
        Args:
            surfaceType (str):                  'line' or 'circle'
            surfaceEquation (list of floats):   [A,B,C] or [h,k,r]
            closestFlag (boolean):              if true, it will return the point closes to the beam's defined x,y,
                                                else return farther point

        Returns:
            list of floats: [x,y] coordinates 
        """

        if surfaceType == 'circle':
            positions = lineCircleIntersection(surfaceEquation, self._beamEquation)

            if closestFlag:
                return min(positions, key=lambda pos: pointPointDistance(pos, self.initialPos) )
            else:
                return max(positions, key=lambda pos: pointPointDistance(pos, self.initialPos) )

        elif surfaceType == 'line':
            return lineLineIntersection(self._beamEquation, surfaceEquation)

    def refractBeam(self, surfaceType, surfaceEquation, closestFlag, RI1, RI2):
        """ refracts light beam

        Args:
            surfaceType (str):                  'line' or 'circle'
            surfaceEquation (list of floats):   [A,B,C] or [h,k,r]
            closestFlag (boolean):              if true, it will return the point closes to the beam's defined x,y,
                                                else return farther point
            RI1 (float):                        refractive index of incident medium
            RI2 (float):                        refractive index of exit medium

        Returns:
            list of floats:                     [A,B,C]
        """

        point = self.propagateBeam(surfaceType, surfaceEquation, closestFlag)
        
        if RI1 > RI2:
            criticalAngle = np.arcsin(RI2 / RI1)
        else:
            criticalAngle = np.pi
        
        if surfaceType == 'circle':
            normalEquation = lineEquationPoints([surfaceEquation[0], surfaceEquation[1]], point)
            
        elif surfaceType == 'line':
            normalEquation = rotateLineAngle(surfaceEquation, np.pi / 2, point)

        self.incidentAngle = lineLineAngle(self.beamEquation, normalEquation)


        if self.incidentAngle > criticalAngle:
            print('attempt to refract failed, passed on to reflect beam')
            return self.reflectBeam(surfaceType, surfaceEquation, closestFlag)
        else:
            self.ExitAngle = np.arcsin(np.sin(self.incidentAngle)*(RI1 / RI2))

        if lineGradient(normalEquation)[0] == 0:
            if self.gradient > 0:
                return rotateLineAngle(normalEquation, self.ExitAngle, point)
            else:
                return rotateLineAngle(normalEquation, -self.ExitAngle, point)
        elif lineGradient(normalEquation)[0] == np.inf:
            if self.gradient > 0:
                return rotateLineAngle(normalEquation, -self.ExitAngle, point)
            else:
                return rotateLineAngle(normalEquation, self.ExitAngle, point)
        elif lineGradient(normalEquation)[0] > 0:
            if self.gradient > 0:
                if self.gradient > lineGradient(normalEquation)[0]:
                    return rotateLineAngle(normalEquation, self.ExitAngle, point)
                else:
                    return rotateLineAngle(normalEquation, -self.ExitAngle, point)
            else:
                return rotateLineAngle(normalEquation, -self.ExitAngle, point)
        elif lineGradient(normalEquation)[0] < 0:
            if self.gradient > 0:
                return rotateLineAngle(normalEquation, self.ExitAngle, point)
            else:
                if self.gradient > lineGradient(normalEquation)[0]:
                    return rotateLineAngle(normalEquation, self.ExitAngle, point)
                else:
                    return rotateLineAngle(normalEquation, -self.ExitAngle, point)
        else:
            raise Exception


    def reflectBeam(self, surfaceType, surfaceEquation, closestFlag):
        point = self.propagateBeam(surfaceType, surfaceEquation, closestFlag)

        if surfaceType == 'cicle':
            normalEquation = lineEquationPoints(surfaceEquation, point)
        elif surfaceType == 'line':
            normalEquation = rotateLineAngle(surfaceEquation, np.pi / 2, point)

        self.incidentAngle = lineLineAngle(self.beamEquation, normalEquation)

        self.ExitAngle = self.incidentAngle
        rotateAngle = 2*self.incidentAngle

        reflectedLine1 = rotateLineAngle(self.beamEquation, rotateAngle, point)
        reflectedLine2 = rotateLineAngle(self.beamEquation, -rotateAngle, point)

        if lineLineAngle(self.beamEquation, normalEquation) > np.radians(45):
            return max([reflectedLine1, reflectedLine2], key=lambda line: lineLineAngle(line, normalEquation))
        else:
            return min([reflectedLine1, reflectedLine2], key=lambda line: lineLineAngle(line, normalEquation))

    def printBeam(self):
        print('')
        for var in vars(self).keys():
            print(f"""{var}: {vars(self)[var]}""")
        print(f'incidentAngle (deg) = {self.incidentAngle / np.pi * 180}')
        print(f'ExitAngle (deg) = {self.ExitAngle / np.pi * 180}')

    @property
    def name(self):
        return self._name

    @property
    def initialPos(self):
        return self._initialPos


    @property
    def gradient(self):
        return self._gradient

    @property
    def intercept(self):
        return self._intercept

    @property
    def beamEquation(self):
        return self._beamEquation