__author__ = 'Anthony G'


class Node:
    # the node number
    number = 0

    # location in 3 coordinates
    x = 0
    y = 0
    z = 0

    # displacement in the 3 dof
    disp_x = 0
    disp_y = 0
    disp_z = 0

    # sensitivity metrics
    sensitivity = 0
    numberOfElements = 0  # the number of elements contributing to this node
    density = 0.0

    def __init__(self, Nodenumber):
        self.number = Nodenumber

    def setLocation(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def setDisplacement(self, dx, dy, dz):
        self.disp_x = dx
        self.disp_y = dy
        self.disp_z = dz
