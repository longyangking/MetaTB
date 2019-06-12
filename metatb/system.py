import numpy as np 

class Point:
    def __init__(self, position, onsite_func):
        self.position = position
        self.onsite_func = onsite_func

    def get_onsite(self, **kwargs):
        onsite = self.onsite_func(**kwargs)
        return onsite

    def get_position(self):
        return self.position

    def set_position(self, position):
        self.position = position

class Link:
    def __init__(self, start_point, end_point, coupling_func):
        self.start_point = start_point
        self.end_point = end_point
        self.coupling_func = coupling_func

    def get_points(self):
        return self.start_point, self.end_point

    def get_coupling(self, **kwargs):
        return self.coupling_func(**kwargs)

class System:
    