import numpy as np
from sympy import Rational
from mpmath import mp


class CR:
    def __init__(self, number):
        self.number = number

    def get_seq(self, n):
        seq = np.zeros(n, dtype=np.longlong)
        number = self.number
        for i in range(n):
            a = mp.floor(number)
            frac = number - a
            seq[i] = a
            if mp.almosteq(number, a):
                number = 0
            else:
                number = mp.fdiv(1, frac)
        return seq

    def get_approximation(self, n):
        oldp = 1
        p = 0
        oldq = 0
        q = 1
        seq = self.get_seq(n)
        for a in seq[1:]:
            if a == 0:
                break

            temp = p
            p = a*p+oldp
            oldp = temp

            temp = q
            q = a*q+oldq
            oldq = temp

        return seq[0] + Rational(p, q)


if __name__ == '__main__':
    mp.dps = 200
    test_cr = CR(3.75)
    print(test_cr.get_seq(10))
    print(test_cr.get_approximation(10))
    pi_cr = CR(mp.pi)
    print(pi_cr.get_seq(100))
    pi_approx = pi_cr.get_approximation(32)
    print(pi_approx)
    print(mp.fdiv(pi_approx.p, pi_approx.q))
    print(mp.pi)