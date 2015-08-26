__author__ = 'Anthony G'

class ContourPoint:
    x = 0
    y = 0
    z = 0

    contourNumber = 0

    def __init__(self, x_local, y_local, z_local,contourNumber_local):
        self.x = x_local
        self.y = y_local
        self.z = z_local
        self.contourNumber = contourNumber_local

