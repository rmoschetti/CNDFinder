from cnd_finder import *
from sage.all import *

IntersectionMatrixK1=matrix([
    [-2,1,0,0,0,0,0,1,1,0,0,0],
    [1,-2,1,0,0,0,0,0,0,0,0,0],
    [0,1,-2,1,0,0,0,0,0,0,0,0],
    [0,0,1,-2,1,0,0,0,0,0,0,0],
    [0,0,0,1,-2,1,0,0,0,0,0,1],
    [0,0,0,0,1,-2,1,0,0,0,0,0],
    [0,0,0,0,0,1,-2,1,0,0,0,0],
    [1,0,0,0,0,0,1,-2,0,0,0,0],
    [1,0,0,0,0,0,0,0,-2,2,0,0],
    [0,0,0,0,0,0,0,0,2,-2,2,0],
    [0,0,0,0,0,0,0,0,0,2,-2,2],
    [0,0,0,0,1,0,0,0,0,0,2,-2]
])

BasisNumK1=[
    [1,0,0,0,0,0,0,0,0,0,0,0], #R1
    [0,1,0,0,0,0,0,0,0,0,0,0], #R2
    [0,0,1,0,0,0,0,0,0,0,0,0], #R3
    [0,0,0,1,0,0,0,0,0,0,0,0], #R4
    [0,0,0,0,1,0,0,0,0,0,0,0], #R5
    [0,0,0,0,0,1,0,0,0,0,0,0], #R6
    [0,0,0,0,0,0,1,0,0,0,0,0], #R7
    [1,1/2,0,1/2,1,1,1,1,1/2,0,0,1/2], #1/2(2R1+R2+R4+2R5+2R6+2R7+2R8+R9+R12)
    [1,1,1,1,1,1/2,0,1/2,1/2,0,0,1/2], #1/2(2R1+2R2+2R3+2R4+2R5+R6+R8+R9+R12)
    [4/2,3/2,2/2,1/2,0,1/2,2/2,3/2,2/2,0,0,0], #1/2(4R1+3R2+2R3+R4+R6+2R7+3R8+2R9)
]

#Expected execution time on an Intel Core i7-7700HQ CPU @ 2.80 Ghz: 4 seconds
FinalResult=CndFinder(IntersectionMatrixK1,BasisNumK1)