imoprt numpy as np
def f(x, y):
    return -(4*x*y+2)/(x**2) - y**2


def yn1(yn, xn, h):
    f1 = f(xn, yn)
    return yn + h*(f1 + f(xn + h, yn + h*f1))/2


def solv(a, y_a, h):
    sol = [y_a]
    for i in range(int((b-a)//h) + 1):
        sol.append(yn1(sol[-1], a + i*h, h))
    return sol


def er(y1, y2):
    return max(list(map(lambda x, y: abs(x-y), y1, y2[::2])))


n = 10
eps = 0.0001
a, b = 1, 2
y_a = -1
h = (b - a)/10
s2 = solv(a, y_a, h)
i = 1
while True:
    s1 = s2
    s2 = solv(a, y_a, h/2)
    if er(s1, s2) < eps:
        break
    h = h/2
    i += 1
o1 = s1[::2**(i-1)]
o2 = s2[::2**i]
print("y" + str(i-1) + " = " + str(o1))
print("y" + str(i) + " = " + str(o2))
print("diff = " + str([o1[i] - o2[i] for i in range(len(o1))]))
