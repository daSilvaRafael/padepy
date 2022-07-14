import sympy as sp
from scipy.linalg import lu_solve


def r(a, b):
    return sp.Rational(a) / sp.Rational(b)


def taylorCoef(func, n):
        """
    Return a list of random ingredients as strings.

    :param kind: Optional "kind" of ingredients.
    :type kind: list[str] or None
    :raise lumache.InvalidKindError: If the kind is invalid.
    :return: The ingredients list.
    :rtype: list[str]
    """
    coef = []
    for coluna in range(0, n + 1):
        derivada = sp.diff(func, x, coluna)
        derivadaEmZero = derivada.subs(x, 0)
        coef.append(derivadaEmZero / sp.Rational(sp.factorial(coluna)))
    return coef


def HankelMatrix(obj, p, q, decPrecision = False):
    A = sp.zeros(q, q)
    try:
        if type(obj) == list:
            elements = obj.copy()
            elements.sort()
            if p + q + 1 <= len(obj):
                for line in range(0, q):
                    for col in range(0, q):
                        if p - q + col + 1 + line < 0:
                            if decPrecision:
                                A[line, col] = 0.0
                            else:
                                A[line, col] = 0
                        else:
                            an = obj[p - q + col + 1 + line]
                            if decPrecision:
                                an_float = sp.N(an, decPrecision)
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
                    # Para todo o n menor que zero n = 0
                    if p - q + col + 1 + line < 0:
                        if decPrecision:
                            A[line, col] = 0.0
                        else:
                            A[line, col] = 0
                    else:
                        derivative = sp.diff(obj, x, p - q + col + line + 1)
                        derivativeAtZero = derivative.subs(x, 0)
                        factorialN = sp.factorial(p - q + col + line + 1)
                        if decPrecision:
                            an_float = sp.N(derivativeAtZero, decPrecision) / sp.N(
                                factorialN, decPrecision
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


def indepVector(obj, p, q, decPrecision = False):
    a = sp.zeros(1, q)
    try:
        if type(obj) == list:
            elements = obj.copy()
            elements.sort()
            if p + q + 1 <= len(obj):
                for row in range(0, q):
                    an = obj[p + row + 1]
                    if decPrecision:
                        an_float = sp.N(an, decPrecision)
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
                if decPrecision:
                    an_float = sp.N(derivativeAtZero, decPrecision) / sp.N(
                        factorialN, decPrecision
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


def denoCoeffForNum(obj, p, q, decPrecision = False):
    cA = sp.zeros(1, p + 1)
    try:
        if type(obj) == list:
            elements = obj.copy()
            elements.sort()
            if p + q + 1 <= len(obj):
                for col in range(0, p + 1):
                    if p - col < 0:
                        if decPrecision:
                            cA[col] = 0.0
                        else:
                            cA[col] = 0
                    else:
                        an = obj[p - col]
                        if decPrecision:
                            an_float = sp.N(an, decPrecision)
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
                    if decPrecision:
                        cA[col] = 0.0
                    else:
                        cA[col] = 0
                else:
                    derivative = sp.diff(obj, x, p - col)
                    derivativeAtZero = derivative.subs(x, 0)
                    factorialN = sp.factorial(p - col)
                    if decPrecision:
                        an_float = sp.N(derivativeAtZero, decPrecision) / sp.N(
                            factorialN, decPrecision
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


def denoCoeff(A, a, q, decPrecision = False):
    b = sp.zeros(1, q + 1)
    if decPrecision:
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


def numCoeff(obj, p, q, bn, decPrecision = False):
    c = sp.zeros(1, p + 1)
    if p < q:
        for col in range(1, p + 2):
            An = sp.Matrix(denoCoeffForNum(obj, p, q, decPrecision)[-col:])
            Bn = sp.Matrix(bn[0:col])
            c[col - 1] = sp.Transpose(An) * Bn
    else:
        col = 1
        while col <= q + 1:
            An = sp.Matrix(denoCoeffForNum(obj, p, q, decPrecision)[-col:])
            Bn = sp.Matrix(bn[0:col])
            c[col - 1] = sp.Transpose(An) * Bn
            col += 1
            j = -1
        while col <= p + 1:
            An = sp.Matrix(denoCoeffForNum(obj, p, q, decPrecision)[-col:j])
            Bn = sp.Matrix(bn[0:])
            c[col - 1] = sp.Transpose(An) * Bn
            col += 1
            j -= 1
    return c


def pade(p, q, obj, decPrecision = 0):
    if type(decPrecision) != int or type(p) != int or type(q) != int:
        print(
            "Numerator and Denominator degree, and decimal precision must be a integer.",
            f"Inputed values: {p, q, decPrecision}",
        )
        return None
    A = HankelMatrix(obj, p, q, decPrecision)
    if f"{type(A)}" != "<class 'sympy.matrices.dense.MutableDenseMatrix'>":
        print(f"Undifined Hankel matrix A. {A}")
        return None
    elif A.det() != 0:
        """Denominator calculation"""
        Dx = sp.zeros(q + 1, 1)
        for row in range(0, q + 1):
            Dx[row] = x ** (row)
            B = sp.transpose(indepVector(obj, p, q, decPrecision))
            bn = denoCoeff(A, B, q, decPrecision)
            if decPrecision:
                B = sp.N(B, decPrecision)
                bn = sp.N(bn, decPrecision)
            Denominator = sp.expand(sp.simplify(bn * Dx))
            """ Numerator calculation """
            Nx = sp.zeros(p + 1, 1)
            for col in range(0, p + 1):
                Nx[col] = x ** (col)
            cn = numCoeff(obj, p, q, bn, decPrecision)
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
