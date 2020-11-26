import math
h = 0.1
t = 0.1
eps = 0.1
N = int(1/h)
a = 0.5

def simpson(y):
    return (h/3)*(y[0] + y[-1] + 2*sum(y[2:-1:2]) + 4*sum(y[1:-1:2]))


def y_f(a):
    y = [0, a]
    for i in range(2, N):
        y.append(h*h*h*(i-1)*math.sqrt(y[i-1])/2 + 2*y[i-1] - y[i-2])
    return y


def fitness(a):
    return simpson(y_f(a)) - 1


def deriv(f, a):
    return (f(a + t) - f(a))/t


def newton(a):
    return a - fitness(a)/deriv(fitness, a)

while abs(fitness(a)) > eps:
    a = newton(a)

answer = y_f(a)
print("answer = " + str(answer))
print("accuracy = " + str(abs(fitness(a))))
