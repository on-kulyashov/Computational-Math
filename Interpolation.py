import numpy as np
from numpy.linalg import inv

x = [0.17453, 0.5236, 0.87267, 1.22173, 1.5708, 1.91986, 2.26893]
y = [0.000003, 0.00018, 0.00227, 0.01770, 0.09688, 0.40481, 1.37878]
n = len(x)


def f(i, j):
    if j - i > 0:
        return (f(i+1, j) - f(i, j-1))/(x[j] - x[i])
    return y[i]


bin_mult = [np.poly1d([1])]
for i in range(n-1):
    bin_mult.append(bin_mult[-1] * np.poly1d([1, -x[i]]))
p = np.poly1d([0])
for i in range(n):
    p += f(0, i)*bin_mult[i]
q = p.deriv()
print(p)

m = [np.array([[x[i]**3, x[i]**2, x[i], 1], [x[i+1]**3, x[i+1]**2, x[i+1], 1], [3*x[i]**2, 2*x[i], 1, 0],
                 [3*x[i+1]**2, 2*x[i+1], 1, 0]]) for i in range(n-1)]
k = [np.array([p(x[i]), p(x[i+1]), q(x[i]), q(x[i+1])]) for i in range(n-1)]
a = [np.matmul(inv(m[i]), k[i].transpose()) for i in range(n-1)]
s = [np.poly1d(a[i]) for i in range(n-1)]
[print(s[i]) for i in range(n-1)]

while True:
    if not n:
        break
    x_0 = float(input())
    print("Common polynomial at this point " + str(p(x_0)))
    i = 0
    while i < n - 1 and x[i] <= x_0:
        i += 1
    assert i >= 1
    print("Spline at this point " + str(s[i-1](x_0)))
