import sympy as sy
import numpy as np


def memoize_doit(fun):
    def doit(self, **hints):
        if "doit" not in self._memoization_cache:
            self._memoization_cache["doit"] = fun(self, **hints)

        return self._memoization_cache["doit"]

    return doit


def memoize_free_symbols(fun):
    def free_symbols(self):
        if "free_symbols" not in self._memoization_cache:
            self._memoization_cache["free_symbols"] = fun(self)

        return self._memoization_cache["free_symbols"]

    return free_symbols


def memoize_xreplace(fun):
    def _xreplace(self, rule):
        hashable_rule = tuple([
            (key, item)
            for key, item in rule.items()
        ])
        cache = self._memoization_cache.setdefault("_xreplace", dict())
        if hashable_rule not in cache:
            cache[hashable_rule] = fun(self, rule)

        return cache[hashable_rule]

    return _xreplace


class MemoizationMatrixExpr(sy.Expr):
    def __new__(cls, *args, **kwargs):
        new_expr = super().__new__(cls, *args, **kwargs)
        new_expr._memoization_cache = dict()
        return new_expr
    
    @memoize_doit
    def doit(self, **hints):
        return super().doit(**hints)

    @property
    @memoize_free_symbols
    def free_symbols(self):
        return super().free_symbols

    @memoize_xreplace
    def _xreplace(self, rule):
        return super()._xreplace(rule)


def inversion_and_det_3x3(m_in):
    """
    Inverts a symbolic 3x3 matrix and calculates the determinant.
    This can be also accomplished with SymPy, the use of this function is for performance reasons only.

    Parameters
    ----------
    m_in : (3,3) SymPy matrix
        Symbolic 3x3 matrix to be inverted

    Returns
    -------
    minv : (3,3) SymPy matrix
        Inverted symbolic matrix
    determinant : SymPy expression
        Determinant of the matrix
    """
    shape = m_in.shape
    m = m_in.reshape(1, 9)

    minv = sy.Matrix(np.zeros((1, 9), dtype=np.object))
    determinant = m[0, 0] * m[0, 4] * m[0, 8] + m[0, 3] * m[0, 7] * m[0, 2] + m[0, 6] * m[0, 1] * m[0, 5] - m[0, 0] * m[
        0, 5] * m[0, 7] - m[0, 2] * m[0, 4] * m[0, 6] - m[0, 1] * m[0, 3] * m[0, 8]
    determinant_inv = 1 / determinant
    minv[0, 0] = (m[0, 4] * m[0, 8] - m[0, 5] * m[0, 7]) * determinant_inv
    minv[0, 1] = (m[0, 2] * m[0, 7] - m[0, 1] * m[0, 8]) * determinant_inv
    minv[0, 2] = (m[0, 1] * m[0, 5] - m[0, 2] * m[0, 4]) * determinant_inv
    minv[0, 3] = (m[0, 5] * m[0, 6] - m[0, 3] * m[0, 8]) * determinant_inv
    minv[0, 4] = (m[0, 0] * m[0, 8] - m[0, 2] * m[0, 6]) * determinant_inv
    minv[0, 5] = (m[0, 2] * m[0, 3] - m[0, 0] * m[0, 5]) * determinant_inv
    minv[0, 6] = (m[0, 3] * m[0, 7] - m[0, 4] * m[0, 6]) * determinant_inv
    minv[0, 7] = (m[0, 1] * m[0, 6] - m[0, 0] * m[0, 7]) * determinant_inv
    minv[0, 8] = (m[0, 0] * m[0, 4] - m[0, 1] * m[0, 3]) * determinant_inv
    return minv.reshape(shape[0], shape[0])


def invert_2x2(matrix):
    """
    Inverts a symbolic 2x2 matrix.

    Parameters
    ----------
    matrix : (2,2) SymPy matrix
        Symbolic 2x2 matrix to be inverted

    Returns
    -------
    inverse : (2,2) SymPy matrix
        Inverted symbolic matrix
    """
    inverse = sy.eye(2)
    determinant = matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0]
    determinant_inv = 1 / determinant
    inverse[0, 0] = matrix[1, 1] * determinant_inv
    inverse[0, 1] = -1 * matrix[0, 1] * determinant_inv
    inverse[1, 0] = -1 * matrix[1, 0] * determinant_inv
    inverse[1, 1] = matrix[0, 0] * determinant_inv
    return inverse


class Inverse(MemoizationMatrixExpr, sy.Inverse):
    __doc__ = sy.Inverse.__doc__

    @memoize_doit
    def doit(self, **hints):
        arg = self.arg  # type: sy.MatrixBase
        arg = arg.doit(**hints)

        if arg.shape == (2, 2):
            inverse = invert_2x2(arg)
            # TODO:
            pass
        elif arg.shape == (3, 3):
            inverse = inversion_and_det_3x3(arg)
            # TODO:
            pass
        else:
            inverse = arg.inv(method="LU")

        return inverse


class Determinant(MemoizationMatrixExpr, sy.Determinant):
    __doc__ = sy.Determinant.__doc__

    @memoize_doit
    def doit(self, **hints):
            arg = self.arg
            if hints.get('deep', True) and isinstance(arg, sy.Basic):
                arg = arg.doit(**hints)
            _eval_determinant = getattr(arg, '_eval_determinant', None)
            if _eval_determinant is not None:
                result = _eval_determinant()
                if result is not None:
                    return result

            return self


class MatMul(MemoizationMatrixExpr, sy.MatMul):
    __doc__ = sy.MatMul.__doc__

    def __new__(cls, *args, **kwargs):
        expanded_args = list()
        for arg in args:
            if isinstance(arg, sy.MatMul):
                expanded_args.extend(arg.args)
            else:
                expanded_args.append(arg)

        return super().__new__(cls, *expanded_args, **kwargs)


class MatAdd(MemoizationMatrixExpr, sy.MatAdd):
    __doc__ = sy.MatAdd.__doc__

    def __new__(cls, *args, **kwargs):
        expanded_args = list()
        for arg in args:
            if isinstance(arg, sy.MatAdd):
                expanded_args.extend(arg.args)
            else:
                expanded_args.append(arg)

        return super().__new__(cls, *expanded_args, **kwargs)


class Transpose(MemoizationMatrixExpr, sy.Transpose):
    __doc__ = sy.Transpose.__doc__
    pass


def main():
    A = sy.MatrixSymbol("A", 2, 2)
    B = Inverse(A.as_explicit())
    C = MatMul(B, B)
    CC = MatMul(C, C)
    D = C.as_explicit()

    C2 = C.xreplace({A[0, 1]: 0, A[1, 0]: 0})
    D2 = C2.as_explicit()
    print("ok")


if __name__ == '__main__':
    main()
