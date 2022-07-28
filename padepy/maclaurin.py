import sympy as sp


def rational_number(a, b):
    """Function to set rational number.

    :param a: numerator.
    :type a: int

    :param b: denominator.
    :type b: int


    :returns: a rational number.
    :rtype: sympy.core.numbers.Rational
    """

    return sp.Rational(a) / sp.Rational(b)


def coefficients(n, var, func):
    """Function to calculate the first n coefficients
    of the Maclaurin expansion.

    :param n: number of coefficients.
    :type n: int

    :param var: function variable.
    :type var: sympy.core.symbol.Symbol

    :param func: User defined function, or sympy bulti-in function like
        sin, cos, log, and exp.
    :type func: sympy.core.mul.Mul or sympy bulti-in function type


    :returns: list with n coefficients.
    :rtype:  list
    """

    coef = []
    for col in range(0, n + 1):
        derivative = sp.diff(func, var, col)
        derivativeAtZero = derivative.subs(var, 0)
        coef.append(derivativeAtZero / sp.Rational(sp.factorial(col)))
    return coef


def polynomial(p, q, var, obj, float_precision=0):
    """Function to calculate the Maclaurin polynomial of order n = p + q

    :param p: degree of the Padé approximant numerator.
    :type p: int

    :param q: degree of the Padé approximant denominator.
    :type q: int

    :param var: function variable.
    :type var: sympy.core.symbol.Symbol

    :param  obj: The obj can be a list of real number, a user defined function,
        or sympy bulti-in function like sin, cos, log, and exp.
    :type obj: list, sympy.core.mul.Mul, or sympy bulti-in function

    :param  float_precision: floating point precision. Default value 0
           for infinite (algebric) precicion.
    :type float_precision: int


    :returns: Maclaurin polynomial of order n = p + q.
    :rtype: sympy.core.add.Add
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
                    fx[col] = var**col
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
                fx[col] = var**col
                derivative = sp.diff(obj, var, col)
                derivativeAtZero = derivative.subs(var, 0)
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


if __name__ == "__main__":
    var = sp.Symbol("x")
    pol = polynomial(3, 3, var, sp.exp(var))
    print()
    print(pol)
    print()
