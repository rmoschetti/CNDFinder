from cnd_finder import *
from sage.all import *

IntersectionMatrixK7=matrix([
     [-2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0], #R1 (E1 in Dolgachev-Kondo)
     [1, -2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2], #R2 (E2 in Dolgachev-Kondo)
     [0, 1, -2, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 2, 0], #R3 (E3 in Dolgachev-Kondo)
     [0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 2, 0, 0, 0], #R4 (E4 in Dolgachev-Kondo)
     [0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 2], #R5 (E5 in Dolgachev-Kondo)
     [0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 2, 0], #R6 (E6 in Dolgachev-Kondo)
     [0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 2, 0, 0], #R7 (E7 in Dolgachev-Kondo)
     [0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 2], #R8 (E8 in Dolgachev-Kondo)
     [1, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 1, 0, 0, 1, 0, 0, 0, 2, 0], #R9 (E9 in Dolgachev-Kondo)
     [0, 1, 1, 0, 0, 0, 0, 0, 0, -2, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0], #R10 (E10 in Dolgachev-Kondo)
     [0, 0, 0, 0, 1, 1, 0, 0, 0, 1, -2, 1, 0, 0, 0, 2, 0, 0, 0, 0], #R11 (E11 in Dolgachev-Kondo)
     [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, -2, 0, 0, 0, 0, 2, 0, 0, 0], #R12 (E12 in Dolgachev-Kondo)
     [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, -2, 0, 0, 2, 0, 0, 0, 0], #R13 (E13 in Dolgachev-Kondo)
     [1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, -2, 0, 0, 2, 0, 0, 0], #R14 (E14 in Dolgachev-Kondo)
     [1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, -2, 0, 0, 2, 0, 0], #R15 (E15 in Dolgachev-Kondo)
     [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, -2, 2, 2, 2, 2], #R16 (K1 in Dolgachev-Kondo)
     [0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 2, -2, 2, 2, 2], #R17 (K2 in Dolgachev-Kondo)
     [0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 2, 2, 2, -2, 2, 2], #R18 (K3 in Dolgachev-Kondo)
     [0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 2, 2, -2, 2], #R19 (K4 in Dolgachev-Kondo)
     [0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, -2]  #R20 (K5 in Dolgachev-Kondo)
])


BasisNumK7=[
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #R1
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #R2
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #R3
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #R4
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #R5
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #R6
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #R7
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #R9
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], #R11
    [1/2, 1/2, 1/2, 1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/2, 0, 0, 0, 0, 0] #1/2(R1+R2+R3+R4+R15)
    ]

#Expected execution time on an Intel Core i7-7700HQ CPU @ 2.80 Ghz: 185 seconds
FinalResult=CndFinder(IntersectionMatrixK7,BasisNumK7)