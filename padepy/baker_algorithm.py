import sympy as sp
from padepy.maclaurin import polynomial


def pade(p, q, var, obj=[], float_precision=0, not_full_path=True):
    """Function to calculate the [p/q](x) Padé approximant using the Baker's recursive algorithm.

    :param p: degree of the Padé approximant numerator.
    :type p: int

    :param q: degree of the Padé approximant denominator.
    :type q: int

    :param var: function variable.
    :type var: sympy.core.symbol.Symbol

    :param  obj: The obj can be a list of real number, a user defined function,
        or sympy bulti-in function like sin, cos, log, and exp. The default empty list
        allows to calculate the Padé approximants generic expressions.
    :type obj: list, sympy.core.mul.Mul, or sympy bulti-in function type

    :param  float_precision: floating point precision. Default value 0 for infinite (algebric) precicion.
    :type float_precision: int

    :param  not_full_path: if set to True the Padé approximants from [p+q/0](x)
    to [0/p+q](x) will be ostructed. If False, from [p+q/0](x) to [p/q](x).
    :type not_full_path: bool


    :returns: The first return is the [p/q](x) Padé approximant. The second is a matrix with all
        Padé approximants contruct until [p/q](x), or all Padé aproximants from [p+q/0](x)
        to [0/p+q](x) if not_full_path set to True.
    :rtype: tuple (sympy.core.mul.Mul, sympy.matrices.dense.MutableDenseMatrix)
    """

    if type(float_precision) != int or type(p) != int or type(q) != int:
        print(
            "Numerator and Denominator degree, and decimal precision must be a integer.",
            f"Inputed values: {p, q, float_precision}",
        )
        return (None, None)
    n = p + q
    if obj == []:
        obj = [f"a{i}" for i in range(0, n + 1)]
    max_iter = 2 * n - 1
    Numerator = sp.zeros(2 + max_iter, 1)
    Denominator = sp.zeros(2 + max_iter, 1)
    f_n = polynomial(p, q, var, obj, float_precision)
    if f_n is None:
        return None
    Numerator[0] = f_n
    if n == 0:
        return (Numerator[0], Numerator[0])

    N0 = sp.Poly(Numerator[0], var)
    if N0.coeff_monomial(var ** (n)) == 0:
        print(f"O(N({var})) < O({var}^{n}) for [{n},0]. Pade's table is not normal.")
        return (None, None)
    Denominator[0] = 1
    f_n_1 = polynomial(p - 1, q, var, obj, float_precision)
    Numerator[1] = f_n_1
    N1 = sp.Poly(Numerator[1], var)
    if N1.coeff_monomial(var ** (n - 1)) == 0:
        print(f"O(N({var})) < O({var}^{n - 1}) for [{n - 1},0]. Pade's table is not normal.")
        return (None, None)
    Denominator[1] = 1
    Pades = sp.zeros(2 + max_iter, 1)
    Pades[0] = sp.Poly(Numerator[0], var) / sp.Poly(Denominator[0], var)
    Pades[1] = sp.Poly(Numerator[1], var) / sp.Poly(Denominator[1], var)
    if p == n and q == 0:
        padePosition = 0
        print(f"Pade [{p},{q}]({var}) stored at matrix index {padePosition}.")
        return (Pades[padePosition], Pades[padePosition])
    i = 2
    j = 2
    m = 1
    padeIndex = False
    while i < j + max_iter:
        if n - m == p and m == q:
            padeIndex = i
        penNumCoeffs = sp.Poly(Numerator[i - 2], var)
        c0 = penNumCoeffs.coeffs()[0]
        lastNumCoeffs = sp.Poly(Numerator[i - 1], var)
        c1 = lastNumCoeffs.coeffs()[0]
        # Baker recursive expressions to construct the approximant [p-m/m]
        Numerator[i] = (1 / c1) * sp.expand(
            sp.simplify((((c1) * Numerator[i - 2]) - (c0) * var * Numerator[i - 1]))
        )
        Denominator[i] = (1 / c1) * sp.expand(
            sp.simplify(((c1) * Denominator[i - 2] - (c0) * var * Denominator[i - 1]))
        )
        Ni = sp.Poly(Numerator[i], var)
        if Ni.coeff_monomial(var ** (n - m)) == 0 and n > m:
            print(
                f"O(N({var})) < O({var}^{n - m}) for [{n - m},{m}]. Pade's table is not normal."
            )
            if padeIndex:
                return (Pades[padeIndex], Pades)
            else:
                return (None, Pades)
        Di = sp.Poly(Denominator[i], var)
        if Di.coeff_monomial(var ** (m)) == 0:
            print(f"O(D({var})) < O({var}^{m}) for [{n - m},{m}]. Pade's table is not normal.")
            if padeIndex:
                return (Pades[padeIndex], Pades)
            else:
                return (None, Pades)
        Pades[i] = sp.Poly(Numerator[i], var) / sp.Poly(Denominator[i], var)
        if padeIndex:
            if not_full_path:
                print(f"Pade [{p},{q}]({var}) stored at matrix index {padeIndex}.")
                pathMatrix = sp.Matrix(Pades[0: padeIndex + 1])
                return (Pades[padeIndex], pathMatrix)
        i += 1
        if i < j + max_iter:
            if n - m - 1 == p and m == q:
                padeIndex = i + 1
            penNumCoeffs = sp.Poly(Numerator[i - 2], var)
            c0 = penNumCoeffs.coeffs()[0]
            lastNumCoeffs = sp.Poly(Numerator[i - 1], var)
            c1 = lastNumCoeffs.coeffs()[0]
            c2 = c1 - c0
            if c2 == 0:
                print(
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
            Ni = sp.Poly(Numerator[i], var)
            if Ni.coeff_monomial(var ** (n - m - 1)) == 0 and n > m + 1:
                print(
                    f"O(N({var})) < O({var}^{n - m - 1}) for [{n - m - 1},{m}]. Pade's table is not normal."
                )
                if padeIndex:
                    return (Pades[padeIndex], Pades)
                else:
                    return (None, Pades)
            Di = sp.Poly(Denominator[i], var)
            if Di.coeff_monomial(var ** (m)) == 0:
                print(
                    f"O(D({var})) < O({var}^{n - m - 1}) for [{n - m - 1},{m}]. Pade's table is not normal."
                )
                if padeIndex:
                    return (Pades[padeIndex], Pades)
                else:
                    return (None, Pades)
            Pades[i] = sp.Poly(Numerator[i], var) / sp.Poly(Denominator[i], var)
            if padeIndex:
                if not_full_path:
                    print(f"Pade [{p},{q}]({var}) stored at matrix index {padeIndex}.")
                    pathMatrix = sp.Matrix(Pades[0: padeIndex + 1])
                    return (Pades[padeIndex], pathMatrix)
        i += 1
        m += 1
    print(f"Pade [{p},{q}]({var}) stored at matrix index {padeIndex}.")
    return (Pades[padeIndex], Pades)


if __name__ == "__main__":
    var = sp.Symbol("x")
    pd, path = pade(3, 3, var, sp.exp(var))
    print()
    print(pd)
    print()
