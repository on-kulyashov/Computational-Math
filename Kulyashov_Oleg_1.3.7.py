import numpy as np

# input
g_0 = 3
r_0 = 21.80089
u_0 = 5.322389*10**5
p_0 = 5.176613*10**12
c_0 = 1.972*10**5
g_3 = 3
r_3 = 10**2
u_3 = 10**6
p_3 = 0
c_3 = 1.972*10**5
# polynomial construction
h_0 = (g_0 + 1) / (g_0 - 1)
h_3 = (g_3 + 1) / (g_3 - 1)
al_0 = r_0 * c_0**2 * (h_0 - 1)
al_3 = r_3 * c_3**2 * (h_3 - 1)
v = u_0 - u_3
x = p_3 / (r_0 * r_3 * v**2)
be_0 = al_0 / (r_0 * r_3 * v**2)
be_3 = al_3 / (r_0 * r_3 * v**2)
de_0 = (h_0 - 1) * r_3
de_3 = (h_3 - 1) * r_0
e = p_0 / (r_0 * r_3 * v**2)
a_6 = (de_0*h_3 - de_3*h_0)**2
a_5 = 2*(h_3 * de_0**2 * (be_3 - h_3 *e) + h_0 * de_3 **2 * (be_0 - x*h_0) - h_0 * h_3 * (de_0*h_3 + de_3*h_0) - de_0*de_3*(h_0*be_3 + h_3*be_0) + de_0*de_3*h_0*h_3*(x+e))
a_4 = (h_0*h_3)**2 + de_0**2 * (h_3**2 * e**2 + be_3**2 - 4*be_3 * h_3 * e) + de_3**2 * (be_0**2 + h_0**2 *x**2 - 4*be_0*x*h_0) - 2*de_0*h_3*(2*be_3*h_0 + h_3*be_0 - e*h_0*h_3) - 2*de_3*h_0*(2*be_0*h_3 + h_0 * be_3 - h_0*h_3*x) - 2*de_0*de_3*(x*e*h_0*h_3 + be_0*de_3) + 2*de_0*de_3*(x+e)*(be_3*h_0+h_3*be_0)
a_3 = 2*(h_0*h_3*(be_0*h_3+be_3*h_0) + de_0**2 * e*be_3*(h_3*e - be_3) + (x*h_0 - be_0)*de_3**2 *be_0*x - de_0*(h_0*be_3**2 - e*be_0*h_3**2) - 2*de_0*h_3*be_3*(be_0 - e*h_0) - de_3*(h_3*be_0**2 - x*be_3*h_0**2) - 2*de_3*be_0*h_0*(be_3 - x*h_3) - e*de_0*de_3*x*(h_0*be_3+h_3*be_0) + de_0*de_3*be_0*be_3*(x+e))
a_2 = be_0**2 * h_3**2 + be_3**2 * h_0**2 + 4*be_0*be_3*h_0*h_3 + de_0**2 * be_3**2 * e**2 + de_3**2 *be_0**2 *x**2 - 2*de_0*(be_0*be_3**2 - e*be_3**2*h_0 - 2*e*h_3*be_0*be_3) - 2*de_3*be_0*(be_0*(be_3 - x*h_3) - 2*be_3*x*h_0) - 2*e*de_0*de_3*be_0*be_3*x
a_1 = 2*(be_0*be_3*(be_0*h_3 + be_3*h_0) + be_0*be_3*(de_0*e*be_3 + de_3*x*be_0))
a_0 = (be_0*be_3)**2
p = np.poly1d([a_6, a_5, a_4, a_3, a_2, a_1, a_0])

# Sturm`s system
def st_syst(p):
    s = [p, np.poly1d.deriv(p)]
    i = s[-2]
    j = s[-1]
    while j.o > 0:
        s.append(-(i / j)[1])
        i = s[-2]
        j = s[-1]
    return s

# number of alterations
def w(c, s):
    k = list(map(lambda q: q(c), s))
    z = 0
    sgn = k[0]
    for i in k:
        if not sgn == np.sign(i):
            sgn = np.sign(i)
            z += 1
    return z

# number of roots on interval
def num_roots(q, a, b):
    s = st_syst(q)
    return w(a, s) - w(b, s)

# bounds for roots
def min_max_r(p):
    m = max(map(abs, p.c[1:]))
    n = max(map(abs, p.c[:-1]))
    max_r = 1 + m / abs(p.c[0])
    min_r = abs(p.c[-1]) / (abs(p.c[-1]) + n)
    return (min_r, max_r)

# intervals with exactly one root
def intervals(p, a, b):
    n = num_roots(p, a, b)
    if n > 1:
        sum_inter = intervals(p, a, (a+b)/2)
        sum_inter.extend(intervals(p, (a+b)/2, b))
        return sum_inter
    elif n:
        return [(a, b)]
    return []

# iterative localization with Sturm`s help
def root(p, interv, precis):
    a, b = interv[0], interv[1]
    leng = b - a
    mid = (a + b) / 2
    if leng > precis:
        if p(a)*p(mid) < 0:
            return root(p, (a, mid), precis)
        return root(p, (mid, b), precis)
    return mid

# all real roots
def roots(p, precision=0.0001):
    rts = []
    for i in intervals(p, *min_max_r(p)):
        rts.append(root(p, i, precision))
    return rts

y = roots(p)
# assert len(y) == 1
'''
Coefficients from the task yield two roots.
One of them leads to nans for shock waves` velocities,
so the other one is taken here.
'''
y = y[1]
# constructing answer
p_1 = y*r_0*r_3*v**2
p_2 = p_1
r_1 = r_0*(p_1*h_0 + al_0)/(p_1 + al_0)
r_2 = r_3*(p_1*h_3 + al_3)/(p_1 + al_3)
u_1 = u_0 + np.sqrt((p_1 - p_0)*(r_1 - r_0)/(r_1*r_0))
u_2 = u_1
d_0 = (r_1*u_1 - r_0*u_0)/(r_1 - r_0)
d_3 = (r_3*u_3 - r_2*u_2)/(r_3 - r_2)
print("d_0 = " + str(d_0) + ", d_3 = " + str(d_3))