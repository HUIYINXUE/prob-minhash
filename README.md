# prob-minhash using f2py

Implement some prob-minhash algorithms via f2py. Algorithm see: O. Ertl, "ProbMinHash - A Class of Locality-Sensitive Hash Algorithms for the (Probability) Jaccard Similarity" in IEEE Transactions on Knowledge & Data Engineering, vol. , no. 01, pp. 1-1, 5555.
doi: 10.1109/TKDE.2020.3021176
keywords: {clustering algorithms;device-to-device communication;estimation error;indexes;measurement;time complexity}
url: https://doi.ieeecomputersociety.org/10.1109/TKDE.2020.3021176

`randgen.f` is a copy from https://www.ucl.ac.uk/~ucakarc/work/randgen.html, but modified to allow call subroutine updating generator seed for more than once.

`my_rand.f` learns from `randgen.f`

Run `f2py -c -m p_minhash p_minhash_0.f90 p_minhash_1.f90 p_minhash_2.f90 p_minhash_3.f90 p_minhash_4.f90 my_rand.f --fcompiler=gnu95 --compiler=unix` to build the module first.