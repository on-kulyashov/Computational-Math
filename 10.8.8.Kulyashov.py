import math


h = 0.01
N = 1000

def y(y1, y0=1):
    yhead = [y0, y1]
    for i in range(N):
        yhead.append(yhead[-2] + 2*h*(-2*yhead[-1] + 1))
    return yhead


yf = list(enumerate(y((1-2*h+math.sqrt(4*(h**2)+1))/2)))
ys = list(enumerate(y((math.exp(-2*h)+1)/2)))
print("With y1 =(1-2h+sqrt(4h^2+1))/2 : " + str(yf))
print("With y1 =(exp(-2*h)+1)/2 : " + str(ys))
