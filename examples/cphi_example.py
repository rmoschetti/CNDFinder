#We compute the combinatorial Phi-invariants of two divisors on
#the Enriques surface with finite automorphism group of type VI.
#In this case, the combinatorial Phi-invariant coincides with the Phi-invariant of the divisor

from cnd_finder import *
from sage.all import *

IntersectionMatrixK6=matrix([
     [-2, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2],
     [1, -2, 1, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 1, -2, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0],
     [0, 0, 1, -2, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0],
     [0, 0, 0, 1, -2, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0],
     [1, 0, 0, 0, 1, -2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0],
     [0, 1, 0, 0, 1, 0, -2, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 1, -2, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 1, 0, 0, 1, 0, 1, -2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0],
     [1, 0, 0, 1, 0, 0, 0, 1, 0, -2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 2, 0, 0, -2, 1, 1, 0, 0, 0, 1, 1, 1, 1],
     [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 1, 1, 0, 1, 1, 0],
     [0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 1, -2, 0, 1, 1, 1, 0, 0, 1],
     [0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, -2, 1, 1, 1, 1, 1, 1],
     [0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 1, 1, -2, 1, 0, 1, 0, 1],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 1, 1, 1, -2, 1, 0, 1, 0],
     [0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, -2, 0, 1, 1],
     [0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, -2, 1, 1],
     [0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, -2, 0],
     [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, -2]
])



BasisNumK6=[
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #R1
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #R2
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #R3
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #R4
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #R5
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #R7
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], #R11
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], #R12
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], #R14
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]  #R17
]

DivisorList=[
    [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0],
    [1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3]
]

#The first divisor D1 in the above list is R2+R12+R3+R17. D1 is the sum of two half-fibers with intersection 1.
#This implies that D1 is nef with positive volume, hence big.
#The program will compute that the Phi-invariant of D1 equals 1.
#Therefore, by Theorem 2.4.14 in Enriques surfaces. I, by Cossec, Dolgachev, Liedtke, D1 has at least one base point.

#The second divisor D2 is the sum of the 10 half-fibers
#[R1, R20],[R2, R12],[R3, R17],[R4, R18],[R5, R13],[R6, R19],[R7, R14],
#[R8, R11],[R9, R15],[R10, R16] divided by 3. This sequence realizes nd(S)=10. It is expected that the Phi-invariant of D2 is 3.
#Due to approximation errors in Python, the code may compute the Phi-invariant to be 2.9999999999999996.
#This can be avoided either by forcing SageMath to work over the rationals, or by computing the Phi-invariant
#of [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], which is 9, and then divide by 3.

#Expected execution time on an Intel Core i7-7700HQ CPU @ 2.80 Ghz: 134 seconds
FinalResult=CPhiFinder(matrix(IntersectionMatrixK6),BasisNumK6,DivisorList)