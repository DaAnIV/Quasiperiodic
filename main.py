import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


def contfrac_to_frac(seq):
    ''' Convert the simple continued fraction in `seq`
        into a fraction, num / den
    '''
    num, den = 1, 0
    for u in reversed(seq):
        num, den = den + num*u, num
    return num, den


class RationalOperator:
    def __init__(self, p, q):
        super().__init__()
        self.alpha = sp.Rational(p, q)

    def _get_frac(self, n):
        na = self.alpha * n
        return na - sp.floor(na)

    def get_v(self, n):
        na = self._get_frac(n)
        return 1 if (1 > na >= 1 - self.alpha) else 0

    def get_sub_matrix(self, n, first_index=0):
        data = np.zeros((n, n), dtype=np.float64)
        for i in range(n):
            if i > 0:
                data[i][i - 1] = 1
            data[i][i] = self.get_v(first_index + i)
            if i < n - 1:
                data[i][i + 1] = 1
        return data


bla = RationalOperator(2, 3)
print(bla.get_sub_matrix(5))
