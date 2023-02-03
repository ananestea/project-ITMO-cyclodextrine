import numpy as np
from math import cos, sin

degree=1  # градус
theta = np.deg2rad(degree) # градус смещени

Xrotation_matrix = np.array([[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]])
Yrotation_matrix = np.array([[cos(theta), 0, sin(theta)], [0, 1, 0], [-sin(theta), 0, cos(theta)]])
Zrotation_matrix = np.array([[cos(theta), -sin(theta), 0], [sin(theta), cos(theta), 0], [0, 0, 1]])
g=np.matrix(Xrotation_matrix)*np.matrix(Yrotation_matrix)*np.matrix(Zrotation_matrix)
print(g)

