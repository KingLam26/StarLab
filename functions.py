# =============== import public modules =============== #
from math import radians
from tkinter import E
import numpy as np


# =============== general functions =============== #
def lineGradient(lineEquation):
    """ returns gradient and y-intercept of a line given by Ax + By + C = 0

    Args:
        lineEquation (list of floats):  [A,B,C] coefficients of the general line equation

    Returns:
        list of floats:                 [gradient, y-intercept] of the line, i.e., y=m*x + c
    """
    A,B,C = lineEquation[0], lineEquation[1], lineEquation[2]
    if B == 0:
        return [np.inf, np.nan]
    else:
        return [-(A/B), -C/B]

def lineAngle(lineEquation):
    """ returns angle the line makes with +ve x-axis

    Args:
        lineEquation (list of floats):  [A,B,C] coefficients of the general line equation

    Returns:
        float:                          angle in radians, between 0 and np.pi
    """
    gradient = lineGradient(lineEquation)[0]

    if gradient == np.inf:
        angle = np.pi / 2
    elif gradient > 0:
        angle = np.arctan(gradient)
    elif gradient < 0:
        angle = np.arctan(gradient) + (2*np.pi)

    if angle > 2*np.pi:
        return angle - 2*np.pi
    else:
        return angle

def lineEquationMC(gradient, intercept):
    """ returns coefficients of the generic line equation given gradient and intercepts of line

    Args:
        m (float):          gradient
        intercept (float):  x-intercept or y-intercept, depending on whether it is a vertical line

    Returns:
        list of floats:     coefficients of the generic line equation, [A,B,C]
                            by default, A = 1, else if A = 0, then B = 1.
    """

    if gradient == np.inf:
        return [1,0,-intercept]
    elif gradient == 0:
        A, B = 0, 1
        return [A,B,-B*intercept]
    else:
        A, B = 1, -1/gradient
        return [A,B, -B*intercept]

def lineEquationPoints(point1, point2):
    """ generates coeffcients of the generic line equation based on coordinates of 2 points on line

    Args:
        point1 / point2 (list of floats):   coodinates of two points on the line, [x1, y1]

    Returns:
        list of floats:                     coefficients of the generic line equation, [A,B,C]
                                            by default, A = 1, else if A = 0, then B = 1.
    """
    x1, y1, x2, y2 = point1[0], point1[1], point2[0], point2[1]

    if x1 == x2:
        return lineEquationMC(np.inf, x1)
    else:
        gradient = (y2 - y1) / (x2 - x1)
        return lineEquationMC( gradient, y1-(gradient*x1) )

def pointPointDistance(point1, point2):
    """ calculates Euclidean distance between two points

    Args:
        point1 / point2 (list of floats):   coodinates of two points on the line, [x1, y1]

    Returns:
        float
    """
    x1, y1, x2, y2 = point1[0], point1[1], point2[0], point2[1]
    
    return np.sqrt( (x2-x1)**2 + (y2-y1)**2 )

def lineLineIntersection(line1, line2):
    """ solves for x-y coordinates of the intersection of two straight lines

    Args:
        line1, line2 (list of floats):  coefficients of the generic line equation, [A1, B1, C1]

    Returns:
        list of floats:                 [x, y] coordinate, intersection of the two lines
    """

    A1, B1, C1, A2, B2, C2 = line1[0], line1[1], line1[2], line2[0], line2[1], line2[2]

    x = (B1*C2 - C1*B2) / (A1*B2 - A2*B1)
    
    if B1 != 0:
        y = (- C1 - A1*x) / B1
    else:
        y = (- C2 - A2*x) / B2

    return [x,y]

def lineLineAngle(line1, line2):
    """ solves for the angle between two straight lines, acute angle

    Args:
        line1, line2 (list of floats):  coefficients of the generic line equation, [A1, B1, C1]

    Returns:
        float:              angle in radians
    """

    m1, m2 = lineGradient(line1)[0], lineGradient(line2)[0]

    if m1 == np.inf:
        angle = np.pi / 2 - np.abs(np.arctan(m2))
    elif m2 == np.inf:
        angle = np.pi / 2 - np.abs(np.arctan(m1))
    else:
        angle = np.arctan(np.abs( (m2 - m1) / (1 + m1*m2)))
    
    if angle > np.pi / 2:
        angle = np.pi - angle

    return angle

def rotatePointAboutPoint(origin, point, angle):
    """ rotates a point counterclockwise by a given angle around a given origin

    Args:
        origin (list of floats):    [x,y]
        point (list of floats):     [x,y]
        angle (float):              radians
    
    Returns:
        list of floats:             [x,y]
    """
    ox, oy = origin[0], origin[1]
    px, py = point[0], point[1]

    qx = ox + np.cos(angle)*(px - ox) - np.sin(angle)*(py - oy)
    qy = oy + np.sin(angle)*(px - ox) + np.cos(angle)*(py - oy)
    return [qx, qy]

def rotateLineAngle(line, alpha, fixedPosition):
    """ rotates straight line by angle alpha, counterclockwise, about fixedPosition

    Args:
        line (list of floats):          coefficients of the generic line equation, [A, B, C]
        alpha (float):                  angle in radians, counterclockwise rotation
        fixedPosition (list of floats): fixed point about which rotation is centered upon, [x,y]

    Returns:
        list of floats:     coefficients of the generic line equation, [A,B,C]
                            by default, A = 1
    """
    
    gradient = lineGradient(line)[0]
    if gradient == np.inf:
        vector = np.array([0, 1])
    else:
        vector = np.array([1, gradient])
    point = vector + np.array(fixedPosition)

    rotatedPoint = rotatePointAboutPoint(fixedPosition, point, alpha)
    
    return lineEquationPoints(rotatedPoint, fixedPosition)

def lineCircleIntersection(circleParameters, lineParameters):
    """ Solves (x - h)^2 + (y - k)^2 = r^2 and Ax + By + C = 0

    Args:
        circleParameters (list of floats):      [h,k,r]
                                                h, k:    coordiante of circle centre
                                                r:       radius of circle
        
        lineParameters (list of floats):        coefficients of the generic line equation, [A, B, C]

    Returns:
        list of floats:                         [x1, y1], [x2, y2] coordinates, intersection of line and circle
    """
    
    h, k, r =  circleParameters[0], circleParameters[1], circleParameters[2]
    A, B, C = lineParameters[0], lineParameters[1], lineParameters[2]

    a = 1 + (A**2 / B**2)
    b = - 2*h + (2*A*C / B**2) + (2*A*k / B)
    d = h**2 + k**2 - r**2 + (2*k*C / B) + ((C**2) / (B**2))
    
    x1 = (-b + np.sqrt(b*b-4*a*d)) / (2*a)
    x2 = (-b - np.sqrt(b*b-4*a*d)) / (2*a)

    y1 = (-C-A*x1) / B
    y2 = (-C-A*x2) / B

    return [x1,y1], [x2,y2]

def coorSysTransform(sysTwoOrigin, sysTwoCoor):
    """ transforms planar coordinates, [x,y] 
        from one orthogonal coordinate system to another,
        assumes axes are mutually parallel

    Args:
        sysTwoOrigin (list of floats):              [x,y] planar coordinates of the second coordinate system origin relative to the first's
        sysOneCoor, sysTwoCoor (list of floats):    [x,y] planar coordinates of the point of interest relative to coordinate system's origin

    Returns:
        (list of floats):                           [x,y] planar coordinates
    """

    return [sysTwoOrigin[0] + sysTwoCoor[0], sysTwoOrigin[1]+sysTwoCoor[1]]


# =============== Functions for plano-convex lens =============== #

def calculatePCLensCT(full_h, radius, edge_t, type):
    """ for a plano-convex lens, calculate its center thickness

    Args:
        full_h (float): full height of the plano-convex lens
        radius (float): radius of the plano-convex lens
        edge_t (float): edge thickness of the plano-convex lens
        type (string):  'concave' or 'convex'

    Returns:
        float: centre thickness of the plano-convex lens
    """

    length = np.sqrt( (full_h / 2)**2 + (edge_t)**2)
    angle = np.pi / 2 + np.arctan(edge_t / (full_h / 2) )
    a = 1
    b = -2*length*np.cos(angle)
    d = length**2 - radius**2

    x1 = (-b + np.sqrt(b*b-4*a*d)) / (2*a)
    x2 = (-b - np.sqrt(b*b-4*a*d)) / (2*a)
    
    tc1, tc2 = radius - x1, radius - x2

    if tc1 > 0:
        tc = tc1
    else:
        tc = tc2
    
    if type == 'convex':
        return tc
    elif type == 'concave':
        return edge_t - np.abs(edge_t - tc)