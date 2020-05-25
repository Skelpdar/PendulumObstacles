import numpy as np

def interpolate(y_0, y_1, y_2):
    a_0 = y_0
    a_1 = (y_1-y_0)/(x_1-x_0)
    a_2 = (y_2-y_1)/((x_2-x_1)*(x_2-x_0)) - (y_1-y_0)/((x_1-x_0)*(x_2-x_0))

    return lambda t: a_0 + a_1*(t-t_0)+a_2*(t-t_0)*(t-t_1)


