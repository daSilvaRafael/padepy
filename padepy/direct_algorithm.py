import sympy as sp
from scipy.linalg import lu_solve


def rational_number(a, b):
    """
    Function to set rational number.
    
    Parameters
    ----------
    a: int 
        Numerator

    b: int
        Denominator
    
    Returns
    -------
    sympy.core.numbers.Rational
        a rational number     
    """
    
    return sp.Rational(a) / sp.Rational(b)


def maclaurin_coefficients(func, n):
    """
    Function to calculate the first n coefficients
    of the Maclaurin expansion. 
    
    Parameters
    ----------
    func: int 
        function to apply the Maclaurin expansion 

    n: int
        number of coefficients
    
    Returns
    -------
    sympy.core.numbers.Rational
        a rational number     
    """

    coef = []
    for coluna in range(0, n + 1):
        derivada = sp.diff(func, x, coluna)
        derivadaEmZero = derivada.subs(x, 0)
        coef.append(derivadaEmZero / sp.Rational(sp.factorial(coluna)))
    return coef


def hankel_matrix(obj, p, q, float_precision=False):
    """
    Function to calculate qxq Hankel matrix.
    
    Parameters
    ----------
    obj: list, sympy.core.mul.Mul or sympy bulti-in function  
        The obj parameter can be a list of real number, a user defined function, or 
        sympy bulti-in function like sin, cos, log, and exp.

    p: int
        degree of the Padé approximant numerator

    q: int
        degree of the Padé approximant denominator

    float_precision: int
        floating point precision. Note that the default value 0 
        is for infinite (algebric) precicion.
        
    Returns
    -------
    sympy.matrices.dense.MutableDenseMatrix  
    """

    A = sp.zeros(q, q)
    try:
        if type(obj) == list:
            elements = obj.copy()
            elements.sort()
            if p + q + 1 <= len(obj):
                for line in range(0, q):
                    for col in range(0, q):
                        if p - q + col + 1 + line < 0:
                            if float_precision:
                                A[line, col] = 0.0
                            else:
                                A[line, col] = 0
                        else:
                            an = obj[p - q + col + 1 + line]
                            if float_precision:
                                an_float = sp.N(an, float_precision)
                                A[line, col] = an_float
                            else:
                                A[line, col] = an
                return A
            else:
                print(
                    f"Error: To construct the Hankel matrix of the Pade [{p},{q}]",
                    "the number of input coefficients needs to be equal to",
                    f"p + q + 1 = {p + q + 1}",
                    f"The input list as {len(obj)} coefficients.",
                )
                return None
        elif type(obj) == tuple:
            print(f"obj = {obj}")
            print(
                f"Error: obj type = {type(obj)}. The obj input must be a function or a list of real numbers."
            )
            return None
        else:
            for line in range(0, q):
                for col in range(0, q):
                    if p - q + col + 1 + line < 0:
                        if float_precision:
                            A[line, col] = 0.0
                        else:
                            A[line, col] = 0
                    else:
                        derivative = sp.diff(obj, x, p - q + col + line + 1)
                        derivativeAtZero = derivative.subs(x, 0)
                        factorialN = sp.factorial(p - q + col + line + 1)
                        if float_precision:
                            an_float = sp.N(derivativeAtZero, float_precision) / sp.N(
                                factorialN, float_precision
                            )
                            A[line, col] = an_float
                        else:
                            A[line, col] = derivativeAtZero / sp.Rational(factorialN)
            return A
    except TypeError:
        print(f"obj = {obj}")
        print(
            f"Error: obj type = {type(obj)}. The obj input must be a function or a list of real numbers."
        )
        return None


def indep_variables(obj, p, q, float_precision = False):
    """
    Function to calculate the independent variables b of
    system Ax=b, where A is the Hankel matrix.
    
    Parameters
    ----------
    obj: list, sympy.core.mul.Mul or sympy bulti-in function  
        The obj parameter can be a list of real number, a user defined function, or 
        sympy bulti-in function like sin, cos, log, and exp.

    p: int
        degree of the Padé approximant numerator

    q: int
        degree of the Padé approximant denominator

    float_precision: int
        floating point precision. Note that the default value 0 
        is for infinite (algebric) precicion.
        
    Returns
    -------
    sympy.matrices.dense.MutableDenseMatrix   
    """
    a = sp.zeros(1, q)
    try:
        if type(obj) == list:
            elements = obj.copy()
            elements.sort()
            if p + q + 1 <= len(obj):
                for row in range(0, q):
                    an = obj[p + row + 1]
                    if float_precision:
                        an_float = sp.N(an, float_precision)
                        a[0, row] = (-1.0) * an_float
                    else:
                        a[0, row] = (-1) * an
                return a
            else:
                print(
                    "Error: The number of input coefficients needs to be equal to p + q + 1",
                    f"= {p + q + 1}. The input list as {len(obj)} coefficients.",
                )
                return None
        elif type(obj) == tuple:
            print(f"obj = {obj}")
            print(
                f"Error: obj type = {type(obj)}. The obj input must be a function or a list of real numbers."
            )
            return None

        else:
            for col in range(0, q):
                derivative = sp.diff(obj, x, p + col + 1)
                derivativeAtZero = derivative.subs(x, 0)
                factorialN = sp.factorial(p + col + 1)
                if float_precision:
                    an_float = sp.N(derivativeAtZero, float_precision) / sp.N(
                        factorialN, float_precision
                    )
                    a[0, col] = (-1.0) * an_float
                else:
                    a[0, col] = (-1) * (derivativeAtZero / sp.Rational(factorialN))
            return a
    except TypeError:
        print(f"obj = {obj}")
        print(
            f"Error: obj type = {type(obj)}. The 'obj' input must be a function or list of real numbers."
        )
        return None


def coefficients_for_numerator(obj, p, q, float_precision = False):
    """
    Function to calculate obj Maclaurin coefficient for
    numerator coefficient calculation.
    
    Parameters
    ----------
    obj: list, sympy.core.mul.Mul or sympy bulti-in function  
        The obj parameter can be a list of real number, a user defined function, or 
        sympy bulti-in function like sin, cos, log, and exp.

    p: int
        degree of the Padé approximant numerator

    q: int
        degree of the Padé approximant denominator

    float_precision: int
        floating point precision. Note that the default value 0 
        is for infinite (algebric) precicion.
        
    Returns
    -------
    sympy.matrices.dense.MutableDenseMatrix 
    """

    cA = sp.zeros(1, p + 1)
    try:
        if type(obj) == list:
            elements = obj.copy()
            elements.sort()
            if p + q + 1 <= len(obj):
                for col in range(0, p + 1):
                    if p - col < 0:
                        if float_precision:
                            cA[col] = 0.0
                        else:
                            cA[col] = 0
                    else:
                        an = obj[p - col]
                        if float_precision:
                            an_float = sp.N(an, float_precision)
                            cA[col] = an_float
                        else:
                            cA[col] = an
                return cA
            else:
                print(
                    "Error: The number of input coefficients needs to be equal to p + q + 1",
                    f"= {p + q + 1}. The input list as {len(obj)} coefficients.",
                )
                return None
        elif type(obj) == tuple:
            print(f"obj = {obj}")
            print(
                f"Error: obj type = {type(obj)}. The obj input must be a function or a list of real numbers."
            )
            return None
        else:
            for col in range(0, p + 1):
                if p - col < 0:
                    if float_precision:
                        cA[col] = 0.0
                    else:
                        cA[col] = 0
                else:
                    derivative = sp.diff(obj, x, p - col)
                    derivativeAtZero = derivative.subs(x, 0)
                    factorialN = sp.factorial(p - col)
                    if float_precision:
                        an_float = sp.N(derivativeAtZero, float_precision) / sp.N(
                            factorialN, float_precision
                        )
                        cA[col] = an_float
                    else:
                        cA[col] = derivativeAtZero / sp.Rational(factorialN)
            return cA
    except TypeError:
        print(f"obj = {obj}")
        print(
            f"Error: obj type = {type(obj)}. The obj input must be a function or a list of real numbers."
        )
        return None


def denominator_coefficients(A, a, q, float_precision = False):
    """
    Function to solve the system Ax=b, where A is the Hankel matrix.
    and x the [p/q](x) Padé approximant denominator coefficients. 
    Parameters
    ----------
    obj: list, sympy.core.mul.Mul or sympy bulti-in function  
        The obj parameter can be a list of real number, a user defined function, or 
        sympy bulti-in function like sin, cos, log, and exp.

    p: int
        degree of the Padé approximant numerator

    q: int
        degree of the Padé approximant denominator

    float_precision: int
        floating point precision. Note that the default value 0 
        is for infinite (algebric) precicion.
        
    Returns
    -------
    sympy.matrices.dense.MutableDenseMatrix    
    """

    b = sp.zeros(1, q + 1)
    if float_precision:
        b[0] = 1.0
    else:
        b[0] = 1
    try:
        Bn = A.LUsolve(a)
        for col in range(1, q + 1):
            b[col] = Bn[-col]
        return b
    except ArithmeticError:
        try:
            Bn = lu_solve(A, a)
            for col in range(1, q + 1):
                b[-col] = Bn[col - 1]
            return b
        except ValueError:
            print("The Ab = c system is ill-posed.")
            return None


def numerator_coefficients(obj, p, q, bn, float_precision = False):
    """
    Function to calculate the [p/q](x) Padé approximant numerator
    coefficients. 
    
    Parameters
    ----------
    obj: list, sympy.core.mul.Mul or sympy bulti-in function  
        The obj parameter can be a list of real number, a user defined function, or 
        sympy bulti-in function like sin, cos, log, and exp.

    p: int
        degree of the Padé approximant numerator

    q: int
        degree of the Padé approximant denominator

    float_precision: int
        floating point precision. Note that the default value 0 
        is for infinite (algebric) precicion.
        
    Returns
    -------
    sympy.matrices.dense.MutableDenseMatrix
        qxq Hankel matrix.     
    """    
    c = sp.zeros(1, p + 1)
    if p < q:
        for col in range(1, p + 2):
            An = sp.Matrix(coefficients_for_numerator(obj, p, q, float_precision)[-col:])
            Bn = sp.Matrix(bn[0:col])
            c[col - 1] = sp.Transpose(An) * Bn
    else:
        col = 1
        while col <= q + 1:
            An = sp.Matrix(coefficients_for_numerator(obj, p, q, float_precision)[-col:])
            Bn = sp.Matrix(bn[0:col])
            c[col - 1] = sp.Transpose(An) * Bn
            col += 1
            j = -1
        while col <= p + 1:
            An = sp.Matrix(coefficients_for_numerator(obj, p, q, float_precision)[-col:j])
            Bn = sp.Matrix(bn[0:])
            c[col - 1] = sp.Transpose(An) * Bn
            col += 1
            j -= 1
    return c


def pade(p, q, obj, float_precision = 0):
    """
    Function to construct the [p/q](x) Padé approximant using the 
    direct algorithm.
    
    Parameters
    ----------
    p: int
        degree of the Padé approximant numerator

    q: int
        degree of the Padé approximant denominator

    obj: list, sympy.core.mul.Mul or sympy bulti-in function  
        The obj parameter can be a list of real number, a user defined function, or 
        sympy bulti-in function like sin, cos, log, and exp. 
        The default empty list allows to calculate some Padé approximants expressions with unknown coefficients. 

    float_precision: int
        floating point precision. Note that the default value 0 
        is for infinite (algebric) precicion.

    Returns
    -------
    sympy.core.mul.Mul
        Return the [p/q](x) Padé approximant.
    """

    if type(float_precision) != int or type(p) != int or type(q) != int:
        print(
            "Numerator and Denominator degree, and decimal precision must be a integer.",
            f"Inputed values: {p, q, float_precision}",
        )
        return None
    A = hankel_matrix(obj, p, q, float_precision)
    if f"{type(A)}" != "<class 'sympy.matrices.dense.MutableDenseMatrix'>":
        print(f"Undifined Hankel matrix A. {A}")
        return None
    elif A.det() != 0:
        """Denominator calculation"""
        Dx = sp.zeros(q + 1, 1)
        for row in range(0, q + 1):
            Dx[row] = x ** (row)
            B = sp.transpose(indep_variables(obj, p, q, float_precision))
            bn = denominator_coefficients(A, B, q, float_precision)
            if float_precision:
                B = sp.N(B, float_precision)
                bn = sp.N(bn, float_precision)
            Denominator = sp.expand(sp.simplify(bn * Dx))
            """ Numerator calculation """
            Nx = sp.zeros(p + 1, 1)
            for col in range(0, p + 1):
                Nx[col] = x ** (col)
            cn = numerator_coefficients(obj, p, q, bn, float_precision)
            Numerator = sp.expand(sp.simplify(cn * Nx))
            """ Pade approximant """
            Pade = sp.Poly(Numerator[0], x) / sp.Poly(Denominator[0], x)
        return Pade
    else:
        print(
            f"The Hankel matrix determinant = {A.det()}. The system Ab = c is impossible or undetermined."
        )
        return None




if __name__ == "__main__":
    x = sp.Symbol("x")
    print()
    print(pade(3, 3, sp.exp(x)))
    print()
else:
    x = sp.Symbol("x")