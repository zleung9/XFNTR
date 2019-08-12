import numpy as np
import matplotlib.pyplot as plt

def laser_gun(origin, alpha, wall=None,celling=None):
    '''
        Define a laser starting from a point with a certain angle,
    until it hits a wall/celling. 
        Input: origin: np.array([x,y])
        Returns the x/y of the hitting point.
    '''
    if wall != None:
        x = wall
        y = origin[1] + np.tan(alpha) * (x - origin[0])
    elif celling != None:
        y = celling
        x = origin[0] + 1/np.tan(alpha) * (celling - origin[1])
    else:
        print "please define a stop"
        return None
    return np.array([x, y])


def central_dist(point, line):
    ''' 
    Calculate the distance from a point to a line.
    Line is defined as: A*x + B*y + C = 0.
    Input line takes in (A,B,C), point takes in [x,y]
    '''
    x, y = point[0], point[1]
    A, B, C = line
    d = np.abs(A*x+B*y+C)/np.sqrt(A**2+B**2)
    return d


def surf_points(x,R):
    '''return an array of surface points given x coordinates
       Input: x, 1D array; curvature in m '''
    x = np.array(x)
    theta = -x / R # sin(theta) = theta, no difference
    y = - R * theta**2/2  # cos(theta) = 1 - theta**2/2, no difference
    surface = np.array([x,y])
    return np.transpose(surface), theta

def hit_surface(point, alpha, curvature, span, sh=0):
    '''
        Given a point and angle, calculate its laser interactig with the surface
        input: points, 1-D array
        output: points, 2-D array
        source: http://mathworld.wolfram.com/Circle-LineIntersection.html. 
    Note that this algorithm applies to circle center at (0,0), sh is the sample height relative NOM position.
        
    '''
    # if point[0]==0.: return np.array([0,0])
    R = curvature 
    cell_left = - span / 2
    cell_right = span / 2
    pt1 = np.array([point[0],point[1]+(R-sh)]) # raise the line up by R.
    pt2 = laser_gun(pt1,alpha,wall=-50) # second point for the line
    
    dx = pt2[0] - pt1[0]
    dy = pt2[1] - pt1[1]
    dr = np.sqrt(dx**2+dy**2)
    D = pt1[0] * pt2[1] - pt2[0] * pt1[1]
    delta = (R*dr)**2-D**2
    if delta < 0:
        hit = np.array([np.inf,np.inf])
    else:
        x1 = (D*dy + dy/abs(dy)*dx*np.sqrt(delta)) / dr**2
        y1 = (-D*dx + abs(dy)*np.sqrt(delta)) / dr**2
#     x2 = (D*dy + sign(dy)*dx*np.sqrt(delta)) / dr**2
#     y2 = (D*dx + abs(dy)*np.sqrt(delta)) / dr**2
        hit = np.array([x1,y1-(R-sh)]) # put the line down by R-sh.
        if hit[0] > cell_right or hit[0] < cell_left:
            hit = np.array([hit[0], np.inf])
    return hit