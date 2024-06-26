import math
from random import random

class Point3:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return f"({self.x},{self.y},{self.z})"

    def dist(self, p):
        return 0.0

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return f"({self.x},{self.y})"


def get_point_in_sphere(r):
    x = (random() * 2 * r) - r
    y = math.sqrt(r * r - x * x)
    y = (random() * 2 * y) - y
    vec = (x * x) + (y * y)
    z = math.sqrt(r * r - vec * vec)
    z = (random() * 2 * z) - z

    return Point3(x, y, z)

# Returns a random point within a circle centered at (0,0)
# with radius r
def get_point_in_circle(r):
    # insert code here
    # use `foo = random()` to generate a random number between [0.0,1.0)
    # use `bar = math.sqrt(y)` to calculate the square root of y
    x = (random() * 2 * r) - r
    y = math.sqrt(r * r - x * x)
    y = (random() * 2 * y) - y
    return Point(x, y)


def fuzz(r):
    for i in range(0, 1_000_000):
        p = get_point_in_circle(r)

        eps = 0.000001
        if (p.x * p.x) + (p.y * p.y) > (r * r) + eps:
            raise Exception(f"{i}: point {p} not in circle")

#p = get_point_in_circle(0.75)
#print(p)
fuzz(1.0)