"""
This interface has two main functions:

* : func:`dense_coefficient_matrix` - computes the dense coefficient matrix of a polynomial system
* : func:`groebner_basis` - computes a Groebner basis with M4GB

"""
from __future__ import absolute_import
import os
import shutil
import subprocess
import sys
from datetime import datetime
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence
from sage.rings.monomials import monomials
from sage.matrix.constructor import matrix
from sage.combinat.integer_vector_weighted import WeightedIntegerVectors
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.misc import walltime

MAX_PRIME = 65521

def check_domains(polyseq, coeff_matrix=False, term_order=None):
    r"""
    Checks whether a suitable base ring
    and term order for the polynomial ring was
    chosen.
    
    INPUT:
    
    - ``polyseq`` -- a polynomial sequence.
    
    - ``coeff_matrix`` -- polynomial sequence is given as matrix, True or False.
     Default: False
    
    - ``term_order`` -- if coefficient matrix is supplied term order must be manually specified.
     Default: None
    
    OUTPUT: raises error if tests fail
    """
    if coeff_matrix:
        field = polyseq.base_ring()
    else:
        ring = polyseq.ring()
        field = ring.base_ring()
        term_order = ring.term_order()
    
    if not field.is_prime_field():
        raise NotImplementedError("base ring must be finite prime field")
    
    if field.characteristic() > MAX_PRIME:
        raise NotImplementedError("maximum prime field size is %s" % MAX_PRIME)
    
    if not (str(term_order) in ["Degree reverse lexicographic term order", "Lexicographic term order", "degrevlex", "lex"]):
        raise NotImplementedError("term order not implemented")

def monomial_from_degree_vector(variables, degree_vector):
    r"""
    Compute a monomial from its degree vector.
    
    INPUT:
    
    - ``variables`` - variables of polynomial ring.
    
    - ``degree_vector`` - degree vector of monomial.
    
    OUTPUT: monomial
    """
    mon = 1
    for var, deg in zip(variables, degree_vector):
        mon *= var**deg
    return mon

def dense_coefficient_matrix(polys):
    r"""
    Compute the dense coefficient matrix
    of a polynomial ideal or sequence.
    Supported coefficient fields are finite prime fields of
    size up to ``MAX_PRIME`` `= 65521 < 2^16`.
    
    INPUT:
    
    - ``polys`` -- an ideal or a polynomial sequence, the generators of an
      ideal., or the dense coefficient matrix of the polynomial sequence.
    
    OUTPUT:
    
    - ``mat`` -- the dense coefficient matrix of the polynomials
    
    Examples:
    This example computes the dense coefficient matrix of a polynomial sequence::
        sage: import sys
        sage: sys.path.insert(0, "PATH_TO_M4GB_SAGE/m4gb_sage/m4gb_sage/")
        sage: import m4gb_sage
        sage: P.<a,b,c> = PolynomialRing(GF(127), order="degrevlex")
        sage: I = sage.rings.ideal.Katsura(P)
        sage: mat = m4gb_sage.dense_coefficient_matrix(I.gens())
        sage: mat
    """
    polyseq = PolynomialSequence(polys)
    check_domains(polyseq)
    ring = polyseq.ring()
    field = ring.base_ring()
    
    max_deg = max([poly.degree() for poly in polyseq])
    degs = []
    for d in range(0, max_deg + 1):
        for deg in WeightedIntegerVectors(d, [1] * len(ring.gens())):
            degs.append(deg)
    monoms = sorted([monomial_from_degree_vector(ring.gens(), deg) for deg in degs], reverse=True)
    
    mat = []
    for poly in polyseq:
        row = []
        for mon in monoms[:-1]:
            try:
                index = poly.monomials().index(mon)
                row.append(poly.coefficients()[index])
            except:
                row.append(field(0))
        mat.append(row)
        row.append(poly([field(0)] * len(ring.gens())))
    mat = matrix(field, mat)
    
    index_non_zero = 0
    zero_col = [field(0)] * len(polyseq)
    for i in range(0, mat.ncols()):
        if mat[:,i].list() != zero_col:
            index_non_zero = i
            break
    mat = mat[:,index_non_zero:]
    
    return mat

def generate_input_file(mat, base_ring, n_vars, term_order):
    r"""
    Generates input file for M4GB in ./tmp folder.
    
    INPUT:
    
    - ``mat`` -- dense coefficient matrix of polynomial sequence.
    
    - ``base_ring`` -- base ring of coefficient matrix.
    
    - ``n_vars`` -- number of variables in polynomial ring.
    
    - ``term_order`` -- term order of polynomial ring.
    
    OUTPUT: file name of input file
    """
    if term_order in ["Degree reverse lexicographic term order", "degrevlex"]:
        term_order = "graded reverse lex order"
    elif term_order in ["Lexicographic term order", "lex"]:
        term_order = "lex order"
    
    dt = datetime.now()
    ts = datetime.timestamp(dt)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    file_name = dir_path + "/../tmp/m4gb_input_" + str(ts) + ".in"
    
    if not os.path.exists(dir_path + "/../tmp"):
        os.mkdir(dir_path + "/../tmp")
    
    with open(file_name, "a") as file:
        lines = []
        lines.append("Galois Field : GF(" + str(base_ring.characteristic()) + ")" + "\n")
        lines.append("Number of variables (n) : " + str(n_vars) + "\n")
        lines.append("Number of equations (m) : " + str(mat.nrows()) + "\n")
        lines.append("Seed : 0" + "\n")
        lines.append("Order : " + term_order + "\n")
        lines.append("\n")
        lines.append("*********************" + "\n")
        file.writelines(lines)
        
        lines = []
        for row in range(0, mat.nrows()):
            line = ""
            for col in range(0, mat.ncols()):
                line += str(mat[row,col]) + " "
            if row < mat.nrows() - 1:
                lines.append(line + ";" + "\n")
            else:
                lines.append(line + ";")
        file.writelines(lines)
    
    return file_name

def read_output_file(file_name_out, ring):
    r"""
    Reads the output file of M4GB and translates
    it back int polynomials of the original polynomial ring
    
    INPUT:
    - ``file_name`` - file name of m4gb output file.
    
    - ``ring`` - polynomial ring of input system to M4GB.
    
    OUTPUT:
    - ``polys`` - output of M4GB in original ring.
    """
    out_variables = ["x" + str(i) for i in range(0, len(ring.gens()))]
    ring_out = PolynomialRing(ring.base_ring(), out_variables, order=ring.term_order())
    with open(file_name_out, "r") as file:
        polys = [line.strip() for line in file]
    polys = [ring(ring_out(poly)) for poly in polys]
    return polys

def create_shell_script(file_name, **kwds):
    r"""
    Creates a shell script in tmp folder for execution of M4GB
    
    INPUT:
    - ``file_name`` - file name of M4GB input file.
    
    - ``threads`` -- integer (default: `1`).
      
    - ``verbosity`` -- integer (default: `3`), display progress info.
    
    OUTPUT:
        file name of M$GB execution script.
    """
    kwds.setdefault('threads', 1)
    kwds.setdefault('verbosity', 3)
    
    dir_path = os.path.dirname(os.path.realpath(__file__))
    m4gb_path = dir_path + "/../m4gb"
    run_file_path = dir_path + "/../tmp/run_m4gb.sh"
    lines = []
    lines.append("#!/bin/bash" + "\n")
    lines.append("cd " + m4gb_path + "\n")
    line = "./solver.sh "
    line += "-t " + str(kwds['threads']) + " "
    line += "-v " + str(kwds['verbosity']) + " "
    line += file_name + " " + file_name.replace(".in", ".out")
    lines.append(line)
    with open(run_file_path, "w") as file:
        file.writelines(lines)
    return run_file_path

def groebner_basis(polys, is_matrix=False, ring=None, **kwds):
    r"""
    Compute a Groebner basis of an ideal using M4GB.
    Supported term orders of the underlying polynomial ring are ``degrevlex``, ``lex``.
    Supported coefficient fields are finite prime fields of
    size up to ``MAX_PRIME`` `= 65521 < 2^16`.
    
    INPUT:
    
    - ``polys`` -- an ideal or a polynomial sequence, the generators of an
      ideal, or its dense coefficient matrix.
    
    - ``is_matrix`` -- is input a dense coefficient matrix, True or False, default: False.
    
    - ``ring`` -- polynomial ring, must only be supplied if polys are given as dense coefficient matrix
    
    - ``threads`` -- integer (default: `1`).
      
    - ``verbosity`` -- integer (default: `3`), display progress info.
    
    OUTPUT: the Groebner basis.
    
    EXAMPLES:
    
    This example computes a Gröbner basis with respect to the degree reverse lexicographic term order::
```     sage: import sys
        sage: sys.path.insert(0, "PATH_TO_M4GB_SAGE/m4gb_sage/m4gb_sage/")
        sage: import m4gb_sage
        sage: P.<a,b,c> = PolynomialRing(GF(127), order="degrevlex")
        sage: I = sage.rings.ideal.Katsura(P)
        sage: gb = m4gb_sage.groebner_basis(I, threads=2. verbosity=3)
        sage: ideal(gb).basis_is_groebner(), ideal(gb) == ideal(I)
        (True, True)
    """
    kwds.setdefault('threads', 1)
    kwds.setdefault('verbosity', 3)
    
    if is_matrix:
        mat = polys
        field = polys.base_ring()
        if ring is None:
            raise NotImplementedError("if polynomial system is given as dense coefficient then polynomial ring must be supplied")
        term_order = ring.term_order()
        check_domains(mat, True, term_order)
        t_mat = -1
    else:
        polyseq = PolynomialSequence(polys)
        check_domains(polyseq)
        ring = polyseq.ring()
        field = ring.base_ring()
        term_order = ring.term_order()
        t_mat = walltime()
        mat = dense_coefficient_matrix(polyseq)
        t_mat = walltime(t_mat)
    gb = None
    
    dir_path = os.path.dirname(os.path.realpath(__file__))
    m4gb_path = dir_path + "/../m4gb"
    cwd = os.getcwd()
    
    file_name = generate_input_file(mat, field, len(ring.gens()), term_order)
    file_name_out = file_name.replace(".in", ".out")
    run_file_path = create_shell_script(file_name, **kwds)
    
    print("M4GB" + "\n")
    t_gb = walltime()
    subprocess.run(["sh", run_file_path], stderr=sys.stderr, stdout=sys.stdout)
    t_gb = walltime(t_gb)
    print("\n")
    
    gb = read_output_file(file_name_out, ring)
    
    tmp_path = dir_path + "/../tmp"
    shutil.rmtree(tmp_path)
    
    if t_mat != -1:
        print("Time needed for coefficient matrix computation: \t" + str(t_mat) + " s.")
    print("Time needed for Gröbner basis computation: \t\t" + str(t_gb) + " s.")
    
    return PolynomialSequence(gb)
