# m4gb_sage
## A SageMath interface for M4GB


This packege is a [SageMath](https://www.sagemath.org/) interface to [M4GB](https://github.com/cr-marcstevens/m4gb) and can be installed as Python package for use with [SageMath](https://www.sagemath.org/).

### Installation
**Requirements**: Linux with a recent version of [SageMath](https://www.sagemath.org/) (tested with SageMath 9.6. installed from Source on Ubuntu 20.04).

First, clone the repository from GitHub and then compile M4GB:
````
git clone https://github.com/SteinerMatthias/m4gb_sage && cd m4gb_sage
./setup.sh
````

### Examples
This example computes a Gröbner basis with respect to the degree reverse lexicographic term order.
```
sage: import sys
sage: sys.path.insert(0, "PATH_TO_M4GB_SAGE/m4gb_sage/")
sage: import m4gb_sage # If it fails, alternatively try: import m4gb_sage.m4gb_sage as m4gb_sage
sage: P.<a,b,c> = PolynomialRing(GF(127), order="degrevlex")
sage: I = sage.rings.ideal.Katsura(P)
sage: gb = m4gb_sage.groebner_basis(I, threads=2. verbosity=3)
sage: ideal(gb).basis_is_groebner(), ideal(gb) == ideal(I)
(True, True)
```

This example computes the dense coefficient matrix of a polynomial sequence.
```
sage: import sys
sage: sys.path.insert(0, "PATH_TO_M4GB_SAGE/m4gb_sage/m4gb_sage/")
sage: import m4gb_sage
sage: P.<a,b,c> = PolynomialRing(GF(127), order="degrevlex")
sage: I = sage.rings.ideal.Katsura(P)
sage: mat = m4gb_sage.dense_coefficient_matrix(I.gens())
sage: mat
[  0   0   0   0   0   0   1   2   2 126]
[  1   0   2   0   0   2 126   0   0   0]
[  0   2   0   0   2   0   0 126   0   0]
```
