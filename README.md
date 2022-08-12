# CNDFinder
SageMath code of the function CND-finder associated with the paper "A computational view on the non-degeneracy invariant for Enriques surfaces", by Riccardo Moschetti, Franco Rota, Luca Schaffler.

For an Enriques surface $S$, the non-degeneracy invariant $nd(S)$ retains information on the elliptic fibrations of $S$ and its polarizations. In the paper related with this code, we introduce a combinatorial version of the non-degeneracy invariant which depends on $S$ together with a configuration of smooth rational curves, and gives a lower bound for $nd(S)$. Here, we provide a SageMath code that computes such combinatorial invariant and we apply it in several examples. 

The file *cnd_finder.py* contains the core SageMath code for the program. 
The folder *examples* contains several applications of the code, toghether with the *result* of the program, for the user's convenience.

We identify a new family of nodal Enriques surfaces satisfying nd(S)=10 which are not general and with infinite automorphism group (*cnd_d16_polarized.py*). We obtain lower bounds on nd(S) for the Enriques surfaces with eight disjoint smooth rational curves studied by Mendes Lopes-Pardini (*cnd_LP1.py*, *cnd_LP2.py*). Finally, we recover Dolgachev and Kondō's computation of the non-degeneracy invariant of the Enriques surfaces with finite automorphism group and provide additional information on the geometry of their elliptic fibrations (*cnd_K1.py*, ..., *cnd_K7.py*).
