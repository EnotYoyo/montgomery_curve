""" This contains three different functions for calculating square roots
modulo a prime.  This first is sqrt1(a,p) which calculates a value for
x where x^2 = a (mod p) where p is a prime and where the variable a is a
quadratic residue in GF(p).  The function sqrt1(a,p) is an implementation of
the Cipolla-Lehmer algorithm.  The second square root function sqrt2(a,p) is
is a different square root algorithm based on exponentiation in GF(p^3) that
works only if p = 5 mod 6 or if p = 4 mod 9 or if p = 7 mod 9.
sqrt1(a,p) is based primarily on the function CL(c,b,p) and sqrt2(a,p) is
based primarily on the function S(d,b,p).  The third square root function
is sqrt4(a,p) which is an improved version of sqrt2(a,p).  sqrt4(a,p) works
for all primes and also avoids having to calculate cube roots and is based on
the function S4(d,b,p).
More information on these functions can be found in the paper
"On Calculating Square Roots in GF(p)" written by David S. Knight
"""

""" This calculates g^e (mod n) for integers e,g, and n """


def exp1(e, g, n):
    t = 1
    sq = g
    e1 = e
    while (e1 != 0):
        if (e1 % 2) == 1:
            t = (sq * t) % n
            e1 = (e1 - 1) // 2
        else:
            e1 = e1 // 2
        sq = (sq * sq) % n
    return (t)


""" This multiplies two polynomials in GF(p^2)
that is it calculates c(x) = a(x)*b(x) <mod q(x)>
where a(x) = a[0]x + a[1] and b(x) = b[0]x + b[1] and
q(x) = x^2 + q[0]x + q[1]
"""


def mult2(a, b, q, n):
    t1 = (a[0] * b[1]) % n
    t2 = (a[1] * b[0]) % n
    t1 = (t1 + t2) % n
    t2 = (a[1] * b[1]) % n
    t3 = ((n - 1) * q[0]) % n
    t4 = ((n - 1) * q[1]) % n
    t5 = (a[0] * b[0]) % n
    t3 = (t5 * t3) % n
    t4 = (t5 * t4) % n
    c = [(t1 + t3) % n, (t2 + t4) % n]
    return (c)


""" This multiplies two polynomials in GF(p^3)
that is it calculates c(x) = a(x)*b(x) <mod q(x)>
where a(x) = a[0]x^2 + a[1]x + a[2] and b(x) = b[0]x^2 + b[1]x + b[2]
and q(x) = x^3 + q[0]x^2 + q[1]x + q[2]
"""


def mult3(a, b, q, n):
    t1 = (a[0] * b[2]) % n
    t2 = (a[1] * b[1]) % n
    t3 = (a[2] * b[0]) % n
    t1 = (t1 + t2) % n
    t1 = (t1 + t3) % n
    t2 = (a[1] * b[2]) % n
    t3 = (a[2] * b[1]) % n
    t2 = (t2 + t3) % n
    t3 = (a[2] * b[2]) % n
    t4 = (a[0] * b[1]) % n
    t5 = (a[1] * b[0]) % n
    t4 = (t4 + t5) % n
    t5 = (a[0] * b[0]) % n
    t6 = ((n - 1) * q[2]) % n
    t6 = (t4 * t6) % n
    t7 = (q[0] * q[2]) % n
    t7 = (t5 * t7) % n
    t6 = (t6 + t7) % n
    t10 = (t6 + t3) % n
    t6 = ((n - 1) * q[1]) % n
    t6 = (t4 * t6) % n
    t7 = (q[0] * q[1]) % n
    t8 = ((n - 1) * q[2]) % n
    t7 = (t7 + t8) % n
    t7 = (t5 * t7) % n
    t6 = (t6 + t7) % n
    t11 = (t6 + t2) % n
    t6 = ((n - 1) * q[0]) % n
    t6 = (t4 * t6) % n
    t7 = (q[0] * q[0]) % n
    t8 = ((n - 1) * q[1]) % n
    t7 = (t7 + t8) % n
    t7 = (t5 * t7) % n
    t6 = (t6 + t7) % n
    t12 = (t6 + t1) % n
    c = [t12, t11, t10]
    return (c)


""" exp2 exponentiates in GF(p^2) by calculating (g(x))^e mod <q(x)> """


def exp2(e, g, q, n):
    t = [0, 1]
    sq = g
    e1 = e
    while (e1 != 0):
        if (e1 % 2) == 1:
            t = mult2(sq, t, q, n)
            e1 = (e1 - 1) // 2
        else:
            e1 = e1 // 2
        sq = mult2(sq, sq, q, n)
    return (t)


""" exp3 exponentiates in GF(p^3) by calculating (g(x))^e mod <q(x)> """


def exp3(e, g, q, n):
    t = [0, 0, 1]
    sq = g
    e1 = e
    while (e1 != 0):
        if (e1 % 2) == 1:
            t = mult3(sq, t, q, n)
            e1 = (e1 - 1) // 2
        else:
            e1 = e1 // 2
        sq = mult3(sq, sq, q, n)
    return (t)


""" This is the main function used in the Cipolla-Lehmer algorithm for
calculating square roots mod a prime p """


def CL(c, b, p):
    t1 = (b * b) % p
    t2 = (4 * c) % p
    t2 = (p - t2) % p
    g = (t1 + t2) % p
    e = (p - 1) // 2
    h = exp1(e, g, p)
    s = 1
    if ((h == 0) or (h == 1)):
        s = 0
    e = e + 1
    t1 = ((p - 1) * b) % p
    t2 = c % p
    q = [t1, t2]
    a = [1, 0]
    t3 = exp2(e, a, q, p)
    t = s * t3[1]
    return (t)


""" This algorithm calculates cube roots mod p assuming p is a prime
that = 5 (mod 6) or p = 4 (mod 9) or p = 7 (mod 9)"""


def cuberoot(g, p):
    t = 0
    p = abs(p)
    g = g % p
    if (p % 6 == 5):
        e = (2 * p - 1) // 3
        t = exp1(e, g, p)
    if (p % 9 == 4) or (p % 9 == 7):
        c = ((p - 1) // 3) % 3
        c = 3 - c
        e = (c * (p - 1) // 3 + 1) // 3
        t = exp1(e, g, p)
        t = 1 * t % p
        t1 = exp1(3, t, p)
        if (t1 != g):
            t = 0
    return (t)


""" This determines if a cubic polynomial q(x) is irreducible in GF(p) """


def ir(q, p):
    t = 1
    a = [0, 1, 0]
    b = exp3(p - 1, a, q, p)
    a = exp3(p + 1, b, q, p)
    b = exp3(p, a, q, p)
    if (b == [0, 0, 1]):
        t = 0
    return (t)


""" This calculates the multiplicative inverse of a in GF(p) """


def inverse(a, p):
    t = exp1(p - 2, a, p)
    return (t)


""" This is the main function used in the GF(p^3) algorithm for
calculating square roots mod a prime p where p = 5 (mod 6) or
p = 4 (mod 9) or p = 7 (mod 9)"""


def S(d, b, p):
    if (p < 4):
        p = 1
    t3 = (p % 6)
    flag = 0
    if (t3 % 6 == 5):
        flag = 1
    if (p % 9 == 4) or (p % 9 == 7):
        flag = 1
    if (flag == 0):
        print("error in S function p should = 5 mod 6 or 4 mod 9 or 7 mod 9")
    t = (b * b) % p
    t = (27 * t) % p
    t = (d + t) % p
    t1 = 4 - (p % 4)
    t1 = (t1 * p + 1) // 4
    t1 = ((p - 1) * t1) % p
    t = (t * t1) % p
    a = cuberoot(t, p)
    q = [0, a, b]
    g = [0, 1, 0]
    s = ir(q, p)
    c = exp3(p, g, q, p)
    t1 = (3 * a) % p
    t2 = inverse(c[0], p)
    t = (t1 * t2) % p
    t = s * t
    if (flag == 0):
        t = 0
    if (flag != 0) and (a == 0):
        t = -1
    return (t)


""" This is the Cipolla-Lehmer algorithm. It calculates a square root
of c mod p assuming c is a quadratic residue mod p where p is a prime > 2"""


def sqrt1(c, p):
    m1 = 50
    t = 0
    c1 = c % p
    for i in range(m1):
        y = CL(c1, ((i + 1) % p), p)
        t1 = (y * y) % p
        if (t1 == c1):
            t = y
            break
    return (t)


""" This is GF(p^3) square root algorithm.  It calculates a square root
of d mod p assuming d is a quadratic residue in GF(p) and p is a prime
where p = 5 (mod 6) or p = 4 (mod 9) or p = 7 (mod 9)"""


def sqrt2(d, p):
    m2 = 50
    t = 0
    d1 = d % p
    if (p % 6 == 5) or (p % 9 == 4) or (p % 9 == 7):
        for i in range(m2):
            y = S(d1, (i + 1) % p, p)
            t1 = (y * y) % p
            if (t1 == d1) and (y != -1):
                t = y
                break
    return (t)


def verify(a, b, p):
    a1 = [0, 0]
    q = [0, a, b]
    t1 = ir(q, p)
    if (t1 == 0):
        print("error in verify function")
        print("the polynomial x^3 + {0}x + {1}".format(a, b))
        print("is not irreducible mod", p)
    if (t1 == 1):
        d = (4 * a * a * a + 27 * b * b) % p
        d = (p - d) % p;
        t = S(d, b, p)
        t2 = inverse(18, p)
        d1 = (t2 * t) % p
        t2 = (p - 1) // 2
        d2 = (t2 * b) % p
        g = [d1, d2]
        e = (p * p - 1) // 3
        q1 = [0, 3]
        a1 = exp2(e, g, q1, p)
    return (a1)


""" Creates a r x c matrix and initializes all entries to 0"""


def init1(r, c):
    s = range(r)
    a = [[0] for n in s]
    for i in range(c - 1):
        for j in range(r):
            a[j].append(0)
    return (a)


""" Creates a r x 1 vector and initializes all entries to 0"""


def init2(r):
    s = range(r)
    a = [0 for n in s]
    return (a)


""" Add two vectors a and b together mod n"""


def add(a, b, n):
    t1 = len(a)
    t2 = len(b)
    t = min(t1, t2)
    s = init2(t)
    for i in range(t):
        s[i] = (a[i] + b[i]) % n
    return (s)


""" Recursive algorithm for calculating a^-1 (mod b) assuming gcd(a,b) = 1"""


def inverse2(a, b):
    if (a == 0 or b == 0):
        return (0)
    if (a == 1):
        return (1)
    else:
        t = (b - ((b * inverse2(b % a, a)) // a))
        return (t)


def multcnvl2(a, b, p):
    t1 = len(a)
    t2 = len(b)
    t = min(t1, t2)
    s = 0
    for i in range(t):
        s1 = (a[i] * b[i]) % p
        s = (s + s1) % p
    return (s)


""" This multiplies two polynomials modulo a cubic polynomial q(x)
that is it calculates c(x) = a(x)*b(x) <mod q(x)>
where a(x) = a[0]x^2 + a[1]x + a[2] and b(x) = b[0]x^2 + b[1]x + b[2]
and q(x) = x^3 + q[0]x^2 + q[1]x + q[2].
Unlike mult3 where the coefficients are integers, in this function the
coefficients are complex integers of the form c1*y^2+c2*y+c3 where y is
a root of the cubic polynomial y^3+q1[0]y^2+q1[1]y+q1[2] = 0 (mod n)
"""


def mult3d(a, b, q, q1, n):
    s = [0, 0, n - 1]

    t1 = mult3(a[0], b[2], q1, n)
    t2 = mult3(a[1], b[1], q1, n)
    t3 = mult3(a[2], b[0], q1, n)
    t1 = add(t1, t2, n)
    t1 = add(t1, t3, n)

    t2 = mult3(a[1], b[2], q1, n)
    t3 = mult3(a[2], b[1], q1, n)
    t2 = add(t2, t3, n)

    t3 = mult3(a[2], b[2], q1, n)
    t4 = mult3(a[0], b[1], q1, n)
    t5 = mult3(a[1], b[0], q1, n)
    t4 = add(t4, t5, n)

    t5 = mult3(a[0], b[0], q1, n)
    t6 = mult3(s, q[2], q1, n)
    t6 = mult3(t4, t6, q1, n)

    t7 = mult3(q[0], q[2], q1, n)
    t7 = mult3(t5, t7, q1, n)
    t6 = add(t6, t7, n)
    t10 = add(t6, t3, n)

    t6 = mult3(s, q[1], q1, n)
    t6 = mult3(t4, t6, q1, n)
    t7 = mult3(q[0], q[1], q1, n)
    t8 = mult3(s, q[2], q1, n)
    t7 = add(t7, t8, n)
    t7 = mult3(t5, t7, q1, n)
    t6 = add(t6, t7, n)
    t11 = add(t6, t2, n)

    t6 = mult3(s, q[0], q1, n)
    t6 = mult3(t4, t6, q1, n)
    t7 = mult3(q[0], q[0], q1, n)
    t8 = mult3(s, q[1], q1, n)
    t7 = add(t7, t8, n)
    t7 = mult3(t5, t7, q1, n)
    t6 = add(t6, t7, n)
    t12 = add(t6, t1, n)

    c = [t12, t11, t10]
    return (c)


""" exp3d is an exponentiation algorithm based on mult3d """


def exp3d(e, g, q, q1, n):
    t = [[0, 0, 0], [0, 0, 0], [0, 0, 1]]
    sq = g
    e1 = e
    while (e1 != 0):
        if (e1 % 2) == 1:
            t = mult3d(sq, t, q, q1, n)
            e1 = (e1 - 1) // 2
        else:
            e1 = e1 // 2
        sq = mult3d(sq, sq, q, q1, n)
    return (t)


""" Primality test """


def is_prime(p):
    p = abs(p)
    if (p == 0):
        p = 1
    m = 50
    y = 0
    b = False
    for i in range(2, m):
        t1 = i % p
        t = exp1((p - 1) // 2, i, p)
        if ((t == 1) or (t == (p - 1)) or (t1 == 0)):
            y = y + 1
    if (y == (m - 2)):
        b = True
    if (p == 1):
        b = False
    return (b)


""" Finds a prime p such that p = 1 (mod c) and p > c*m
If c is odd then the prime p = 3 (mod 4)
If 2^a divides c and 2^(a+1) does not divide c where a>0
then p = (2^a)+1 mod (2^(a+1)) """


def findprime(c, m):
    if (c < 1):
        c = 1
    if (m < 1):
        m = 1
    if (c % 2 == 1):
        c = 2 * c
        m = (m + (m % 2)) // 2
    t = False
    i = m + (m + 1) % 2
    while (t == False):
        t1 = c * i + 1
        t = is_prime(t1)
        i = i + 2
    return (t1)


""" This is the main function used in the improved GF(p^3) algorithm for
calculating the square root of a quadratic residue d mod a prime p where p > 3
based on a random variable b were 0 < b < p"""


def S4(d, b, p):
    p = abs(p)
    if (p == 0):
        p = 2
    if (d < 0):
        d = abs(d)
        d = (p - d) % p
    if (b < 0):
        b = abs(b)
        b = (p - b) % p
    d = d % p
    b = b % p
    t = (b * b) % p
    t = (27 * t) % p
    t = (d + t) % p
    t1 = 4 - (p % 4)
    t1 = (t1 * p + 1) // 4
    t1 = ((p - 1) * t1) % p
    t = (t * t1) % p
    q1 = [0, 0, p - t]
    q = [[0, 0, 0], [0, 1, 0], [0, 0, b]]
    g = [[0, 0, 0], [0, 0, 1], [0, 0, 0]]
    t = 0
    t1 = (q1[2] * b) % p

    if (p % 6 == 1) and (t1 != 0):
        s1 = exp3d(p, g, q, q1, p)
        t1 = s1[0][1]
        t = inverse(t1, p)
        t = (3 * t) % p
        t1 = t ** 2 % p
        if (t1 != d):
            t2 = inverse(d, p)
            t2 = (t1 * t2) % p
            t = (t * t2) % p

    if (p % 6 == 5) and (t1 != 0):
        s1 = exp3d(p ** 2, g, q, q1, p)
        t1 = s1[0][1]
        t = inverse(t1, p)
        t = (3 * t) % p
        t = (p - t) % p

    return (t)


""" This is the improved GF(p^3) square root algorithm and is based on the
function S4(d,b,p). It calculates a square root of d mod p assuming d is a
quadratic residue in GF(p).  It returns the same output as sqrt2
if p = 5 (mod 6) or if p = 4 (mod 9) or p = 7 (mod 9) where p is a prime.
But unlike sqrt2, sqrt4 also works in all cases if p = 1 (mod 6)
((sqrt2 doesn't work if p = 1 (mod 9)).  Also sqrt4 improves the efficiency
if p = 1 (mod 6) by avoiding having to calculate a cube root in GF(p)"""


def sqrt4(d, p):
    m2 = 50
    t = 0
    d1 = d % p
    t = is_prime(p)
    if (t == True):
        for i in range(m2):
            y = S4(d1, (i + 1) % p, p)
            t1 = (y * y) % p
            if (t1 == d1):
                t = y
                break
    if (t == False):
        t = 0
    return (t)
