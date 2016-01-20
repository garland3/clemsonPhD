__author__ = 'Anthony G'

class Element:
    def __init__(self, number, body, densityDesignVar, typeNumber, nodesindex):
        self.nodeNumber = number
        self.bodyNumber = body
        self.elementTypeNumber = typeNumber
        self.nodesIndex = nodesindex  # list with nodes in this element
        self.numberNodes = 0
        self.massMatrix = []
        self.stiffnessMatrix = []
        self.mass = 0
        self.densityDesignVar = densityDesignVar
        self.sensitivity = 0