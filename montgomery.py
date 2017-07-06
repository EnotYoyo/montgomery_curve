import gmpy2
import random
from cipolla_lehmer import sqrt1
import time
from sympy import *

# P = 115792089237316193816632940749697632357179977837030573159940048803231482807789
# A = 25145277331988958011270560234050020930586280001375125240020481842714790931915
# B = 69211765927878610661416544348746032225986270934521153725036655380879831582035
# X = 2978887011322666093963289506157390814196814788524229414376294405714137817842
# Y = 20970510501817807007505964668159123466115295583918396465274387362709663676660
# P = 2 ** 255 - 19
P = 70376410737831876759654067040077648027271263871227689830927128611623145214459
A = 486662
B = 1


class Point(object):
    def __init__(self, x=0, y=None):
        self.z = 1
        if x == 0:
            self.x, self.y = 0, 0
        elif x is None:
            while True:
                self.x = random.randint(0, P)
                if gmpy2.jacobi((self.x ** 3 + A * (self.x ** 2) + self.x) * gmpy2.invert(B, P), P) == 1:
                    self.y = int(sqrt1((self.x ** 3 + A * (self.x ** 2) + self.x) * gmpy2.invert(B, P), P)) % P
                    break
        else:
            self.x = x % P
            if gmpy2.jacobi((self.x ** 3 + A * (self.x ** 2) + self.x) * gmpy2.invert(B, P), P) != 1:
                print("Bad x coordinates:", self.x)
                exit(0)
            if y is None:
                self.y = None
            else:
                self.y = y
                tmp1 = (B * self.y ** 2) % P
                tmp2 = (self.x ** 3 + A * (self.x ** 2) + self.x) % P
                if tmp1 != tmp2:
                    print("Bad y coordinates:", y)
                    exit(0)

    def add(self, other, difference=None):
        if self.x == 0 and self.y == 0:
            return other
        if other.x == 0 and other.y == 0:
            return self

        if self.x == other.x:
            if self.y == other.y:
                xx1 = self.x ** 2
                x3 = pow(xx1 - self.z ** 2, 2, P)
                z3 = 4 * self.x * self.z * (xx1 + A * self.x * self.z + self.z ** 2) % P
                z3 = gmpy2.invert(z3, P)
                return Point(int(x3 * z3))

        if difference is not None:
            f1 = (other.x - other.z) * (self.x + self.z) % P
            f2 = (other.x + other.z) * (self.x - self.z) % P
            x3 = difference.z * pow(f1 + f2, 2, P)
            z3 = difference.x * pow(f1 - f2, 2, P)
            return Point(int(x3 * gmpy2.invert(z3, P)))

        if self.y is None:
            self._recovery_y()

        yy1 = (other.y - self.y) % P
        xx1 = (other.x - self.x) % P
        x3 = B * yy1 ** 2 * gmpy2.invert(xx1 ** 2, P) - A - self.x - other.x
        y3 = (2 * self.x + other.x + A) * yy1 * gmpy2.invert(xx1, P) - B * yy1 ** 3 * gmpy2.invert(xx1 ** 3, P) - self.y

        return Point(x3 % P, y3 % P)

    def _recovery_y(self):
        self.y = int(sqrt1((self.x ** 3 + A * (self.x ** 2) + self.x) * gmpy2.invert(B, P), P)) % P
        self.y = min(P - self.y, self.y)

    def mul_binary(self, other):
        r = Point()
        m2 = self
        n = other
        for i in [int(digit) for digit in bin(n)[2:]]:
            if i == 1:
                r = m2.add(r)
            n, m2 = n >> 1, m2.add(m2)
        return r

    def mul(self, other):
        r0 = self
        r1 = self.add(self)
        n = other
        for i in [int(digit) for digit in bin(n)[3:]]:
            if i == 0:
                r1 = r0.add(r1, self)
                r0 = r0.add(r0)
            else:
                r0 = r0.add(r1, self)
                r1 = r1.add(r1)
        return r0

    def __str__(self):
        if self.y is None:
            self._recovery_y()

        return '({}, {})\n({}, {})\n'.format(self.x * gmpy2.invert(self.z, P), self.y * gmpy2.invert(self.z, P),
                                             self.x * gmpy2.invert(self.z, P), P - self.y * gmpy2.invert(self.z, P))


def generate_montgomery_curve():
    while True:
        b = 1  # random.randint(2, P)
        if gmpy2.jacobi(b, P) != 1:
            continue

        x = random.randint(2, P)
        a = random.randint(2, P)
        if b * (a ** 2 - 4) == 0 \
                or a == 2 \
                or a == P - 2 \
                or gmpy2.jacobi((x ** 3 + a * (x ** 2) + x) * gmpy2.invert(b, P), P) != 1:
            continue
        y = sqrt1((x ** 3 + a * (x ** 2) + x) * gmpy2.invert(b, P), P)
        return a, b, x, int(y)


def generate_weierstrass_curve_from_montgomery(ma, mb, mx, my, p):
    wx = (mx * gmpy2.invert(mb, p) + ma * gmpy2.invert(3 * mb, p)) % p
    wy = (my * gmpy2.invert(mb, p)) % p
    wa = ((3 - ma ** 2) * gmpy2.invert(3 * (mb ** 2), p)) % p
    wb = ((2 * (ma ** 3) - 9 * ma) * gmpy2.invert(27 * mb, p)) % p
    assert 0 < wa < p and 0 < wb < p and p > 2
    assert (4 * (wa ** 3) + 27 * (wb ** 2)) % p != 0
    return wa, wb, wx, wy


def generate_montgomery_curve_from_weierstrass(wa, wb, wx, wy, p):
    from sympy.polys.domains import ZZ
    import sympy
    from sympy.polys.galoistools import gf_factor
    f = [1, 0, wa, wb]
    sympy.Poly.from_list(f, sympy.Symbol('x'))
    factor = gf_factor(f, p, ZZ)
    x = None
    for r in factor:  # fucking sympy
        if isinstance(r, list):
            for root in r:
                if len(root) == 2 \
                        and len(root[0]) == 2 \
                        and gmpy2.jacobi(3 * ((p - root[0][-1]) ** 2) + wa, p) == 1:
                    x = p - root[0][-1]
                    break

    assert x is not None

    s = sqrt1(gmpy2.invert(3 * x ** 2 + wa, p), p)
    mb = s
    ma = (3 * x * s) % p
    mx = (s * (wx - x)) % p
    my = (s * wy) % p
    return (ma, p - ma), (mb, p - mb), (mx, p - mx), (my, p - my)


def timer(f):
    def tmp(*args, **kwargs):
        t = time.time()
        n = 10
        for x in range(n):
            res = f(*args, **kwargs)
        print('{name} {time:f}'.format(name=args[1], time=(time.time() - t) / n))
        return res

    return tmp


@timer
def mul(point, length):
    """
        Testing the speed of Montgomery's ladder
        :param point:
        :param length:
        :return:
    """
    return point.mul(random.getrandbits(length))


@timer
def mul_bynary(point, length):
    """
       Testing the speed of binary multiplication
       :param point:
       :param length:
       :return:
    """
    return point.mul_binary(random.getrandbits(length))


def test_arithmetic_operations():
    ma, mb, mx, my = generate_montgomery_curve()
    print('montgomery params:\n ma = {}\n mb = {}\n mx = {}\n my = {}\n'.format(ma, mb, mx, my))
    mul = random.randint(2, 100)
    P0 = Point(mx, my)
    P1 = P0.mul(mul)
    P3 = P0.mul(9010)
    print('{} = {} * {} = '.format(P1, mul, P0))
    print('P1 + 9010 * P0 = ', P1.add(P3, P0))


def test_translation():
    ma, mb, mx, my = generate_montgomery_curve()
    wa, wb, wx, wy = generate_weierstrass_curve_from_montgomery(ma, mb, mx, my, P)

    # print('weierstrass params:\np = {}\nA = {}\nB = {}\nx_G = {}\ny_G = {}\n'.format(P, wa, wb, wx, wy))
    # print('montgomery params:\np = {}\nA = {}\nB = {}\nx_G = {}\ny_G = {}\n'.format(P, ma, mb, mx, my))

    tests = generate_montgomery_curve_from_weierstrass(wa, wb, wx, wy, P)
    assert ma == tests[0][0] or ma == tests[0][1]
    assert mb == tests[1][0] or mb == tests[1][1]
    assert mx == tests[2][0] or mx == tests[2][1]
    assert my == tests[3][0] or my == tests[3][1]


if __name__ == '__main__':
    test_translation()
