from e import interpolate_vec
from e import interpolate
from scipy.optimize import ridder 
import numpy as np

"""
The root is found in the interval t_i-1 to t_i with Ridders method.
Although for a polynomial of degree two, we can find the zeroes analytically.
But here we implement a more general numerical method for posterity.
"""


def has_solution(interpolated_vec,t_last, t_current, alpha_obst):
    try:
        ridder(lambda t: interpolated_vec(t)[0][0]-alpha_obst, t_last, t_current)
        return True
    except:
        return False

