import sympy as sp


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
    of func Maclaurin expansion. 
    
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
    for col in range(0, n + 1):
        derivative = sp.diff(func, x, col)
        derivativeAtZero = derivative.subs(x, 0)
        coef.append(derivativeAtZero / sp.Rational(sp.factorial(col)))
    return coef


def trunc_series(obj, p, q, float_precision=0):
    """
    Function to construct Maclaurin expansion of order n = p + q
    
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
    sympy.core.add.Add
        Maclaurin expansion of order n = p + q     
    """

    n = p + q
    an = sp.zeros(n + 1, 1)
    fx = sp.zeros(1, n + 1)
    try:
        if type(obj) is list:
            elements = obj.copy()
            elements.sort()
            if len(obj) >= n + 1:
                for col in range(0, n + 1):
                    coef = obj[col]
                    fx[col] = x**col
                    if float_precision:
                        an[col] = sp.N(coef, float_precision)
                    else:
                        an[col] = coef
            else:
                print(
                    f"Error: To construct the Pade [{p},{p}]",
                    "the number of input coefficients needs to be equal to",
                    f"p + q + 1 = {p + q + 1}",
                    f"The input list as {len(obj)} coefficients.",
                )
                return None
            Sn = fx * an
            return Sn[0]
        elif type(obj) == tuple:
            print(f"obj = {obj}")
            print(
                f"Error: obj type = {type(obj)}. The obj input must be a function or a list of real numbers."
            )
            return None
        else:
            for col in range(0, n + 1):
                fx[col] = x**col
                derivative = sp.diff(obj, x, col)
                derivativeAtZero = derivative.subs(x, 0)
                factorialN = sp.factorial(col)
                if float_precision:
                    an_float = sp.N(derivativeAtZero, float_precision) / sp.N(
                        factorialN, float_precision
                    )
                    an[col] = an_float
                else:
                    an[col] = derivativeAtZero / sp.Rational(factorialN)
            Sn = fx * an
            return Sn[0]
    except TypeError:
        print(f"obj = {obj}")
        print(
            f"The obj type = {type(obj)}. The obj input must be a function or list of real numbers.",
        )
        return None


def pade(p, q, obj=[], float_precision=0, not_full_path=True):
    """
    Function to construct the [p/q](x) Padé approximant using the 
    Baker's recursive algorithm.
    
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

    not_full_path: bool
        if set to True the Padé approximants from [p/q](x) to [0/p+q](x) will be
        costructed. Otherwise the algorithm stops once the [p/q](x) is constructed.

    Returns
    -------
    tuple(sympy.core.mul.Mul, sympy.matrices.dense.MutableDenseMatrix)
        The first return is the [p/q](x) Padé approximant and
        the second is a matrix with all Padé approcimant contruct until [p/q](x) or
        all Padé aproximants from [p+q/1](x) to [0/p+q](x) if not_full_pathset to
        True.    
    """

    if type(float_precision) != int or type(p) != int or type(q) != int:
        print(
            "Numerator and Denominator degree, and decimal precision must be a integer.",
            f"Inputed values: {p, q, float_precision}",
        )
        return (None, None)
    n = p + q
    if obj == []:
        for i in range(0, n + 1):
            obj.append(f"a{i}")
    max_iter = 2 * n - 1
    Numerator = sp.zeros(2 + max_iter, 1)
    Denominator = sp.zeros(2 + max_iter, 1)
    f_n = trunc_series(obj, p, q, float_precision)
    if f_n is None:
        return None
    Numerator[0] = f_n
    if n == 0:
        return (Numerator[0], Numerator[0])

    N0 = sp.Poly(Numerator[0], x)
    if N0.coeff_monomial(x ** (n)) == 0:
        print(f"O(N(x)) < O(x^{n}) for [{n},0]. Pade's table is not normal.")
        return (None, None)
    Denominator[0] = 1
    f_n_1 = trunc_series(obj, p - 1, q, float_precision)
    Numerator[1] = f_n_1
    N1 = Ni = sp.Poly(Numerator[1], x)
    if N1.coeff_monomial(x ** (n - 1)) == 0:
        print(f"O(N(x)) < O(x^{n - 1}) for [{n - 1},0]. Pade's table is not normal.")
        return (None, None)
    Denominator[1] = 1
    Pades = sp.zeros(2 + max_iter, 1)
    Pades[0] = sp.Poly(Numerator[0], x) / sp.Poly(Denominator[0], x)
    Pades[1] = sp.Poly(Numerator[1], x) / sp.Poly(Denominator[1], x)
    if p == n and q == 0:
        padePosition = 0
        print(f"Pade [{p},{q}](x) stored at matrix index {padePosition}.")
        return (Pades[padePosition], Pades[padePosition])
    i = 2
    j = 2
    m = 1
    padeIndex = False
    while i < j + max_iter:
        if n - m == p and m == q:
            padeIndex = i
        penNumCoeffs = sp.Poly(Numerator[i - 2], x)
        c0 = penNumCoeffs.coeffs()[0]
        lastNumCoeffs = sp.Poly(Numerator[i - 1], x)
        c1 = lastNumCoeffs.coeffs()[0]
        # Baker recursive expressions to construct the approximant [p-m/m]
        Numerator[i] = (1 / c1) * sp.expand(
            sp.simplify((((c1) * Numerator[i - 2]) - (c0) * x * Numerator[i - 1]))
        )
        Denominator[i] = (1 / c1) * sp.expand(
            sp.simplify(((c1) * Denominator[i - 2] - (c0) * x * Denominator[i - 1]))
        )
        Ni = sp.Poly(Numerator[i], x)
        if Ni.coeff_monomial(x ** (n - m)) == 0 and n > m:
            print(
                f"O(N(x)) < O(x^{n - m}) for [{n - m},{m}]. Pade's table is not normal."
            )
            if padeIndex:
                return (Pades[padeIndex], Pades)
            else:
                return (None, Pades)
        Di = sp.Poly(Denominator[i], x)
        if Di.coeff_monomial(x ** (m)) == 0:
            print(f"O(D(x)) < O(x^{m}) for [{n - m},{m}]. Pade's table is not normal.")
            if padeIndex:
                return (Pades[padeIndex], Pades)
            else:
                return (None, Pades)
        Pades[i] = sp.Poly(Numerator[i], x) / sp.Poly(Denominator[i], x)
        if padeIndex:
            if not_full_path:
                print(f"Pade [{p},{q}](x) stored at matrix index {padeIndex}.")
                pathMatrix = sp.Matrix(Pades[0: padeIndex + 1])
                return (Pades[padeIndex], pathMatrix)
        i += 1
        if i < j + max_iter:
            if n - m - 1 == p and m == q:
                padeIndex = i + 1
            penNumCoeffs = sp.Poly(Numerator[i - 2], x)
            c0 = penNumCoeffs.coeffs()[0]
            lastNumCoeffs = sp.Poly(Numerator[i - 1], x)
            c1 = lastNumCoeffs.coeffs()[0]
            c2 = c1 - c0
            if c2 == 0:
                print (
                    f"Baker's algorithm fail to construct [{n - m - 1},{m}].",
                    f"The normalization coefficient with {float_precision} decimal precison is 0 = {c1} - {c0} ,",
                    f"which implies [{n - m - 1},{m}] = inf / inf.",
                )
                return (None, None)
            # Baker recursive expressions to construct the approximant [p-m-1/m]
            Numerator[i] = (1 / c2) * sp.expand(
                sp.simplify(((c1) * Numerator[i - 2]) - (c0) * Numerator[i - 1])
            )
            Denominator[i] = (1 / c2) * sp.expand(
                sp.simplify(((c1) * Denominator[i - 2] - (c0) * Denominator[i - 1]))
            )
            Ni = sp.Poly(Numerator[i], x)
            if Ni.coeff_monomial(x ** (n - m - 1)) == 0 and n > m + 1:
                print(
                    f"O(N(x)) < O(x^{n - m - 1}) for [{n - m - 1},{m}]. Pade's table is not normal."
                )
                if padeIndex:
                    return (Pades[padeIndex], Pades)
                else:
                    return (None, Pades)
            Di = sp.Poly(Denominator[i], x)
            if Di.coeff_monomial(x ** (m)) == 0:
                print(
                    f"O(D(x)) < O(x^{n - m - 1}) for [{n - m - 1},{m}]. Pade's table is not normal."
                )
                if padeIndex:
                    return (Pades[padeIndex], Pades)
                else:
                    return (None, Pades)
            Pades[i] = sp.Poly(Numerator[i], x) / sp.Poly(Denominator[i], x)
            if padeIndex:
                if not_full_path:
                    print(f"Pade [{p},{q}](x) stored at matrix index {padeIndex}.")
                    pathMatrix = sp.Matrix(Pades[0: padeIndex + 1])
                    return (Pades[padeIndex], pathMatrix)
        i += 1
        m += 1
    print(f"Pade [{p},{q}](x) stored at matrix index {padeIndex}.")
    return (Pades[padeIndex], Pades)


if __name__ == "__main__":
    x = sp.Symbol("x")
    print()
    pd, path = pade(3, 3, sp.exp(x))
    print()
    print(pd)
    print()
else:
    x = sp.Symbol("x")