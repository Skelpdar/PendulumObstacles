import numpy as np

def interpolate(y_0, y_1, y_2,x_0,x_1,x_2):
    a_0 = y_0
    a_1 = (y_1-y_0)/(x_1-x_0)
    a_2 = (y_2-y_1)/((x_2-x_1)*(x_2-x_0)) - (y_1-y_0)/((x_1-x_0)*(x_2-x_0))

    return lambda t: a_0 + a_1*(t-x_0)+a_2*(t-x_0)*(t-x_1)

def interpolate_vec(y_0, y_1, y_2, x_0, x_1, x_2):
    return lambda t: np.array([[interpolate(y_0[0][0],y_1[0][0], y_2[0][0], x_0, x_1, x_2)(t)],[interpolate(y_0[1][0],y_1[1][0],y_2[1][0], x_0, x_1, x_2)(t)]])

