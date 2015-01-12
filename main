__author__ = 'hose314'
# -*- coding: utf-8 -*-
import numpy as np
import math
import pandas as pd
from parser import expr
from copy import copy

n = 20
m = 400
l = 5
T = 1
#предпологаяется что функция непрервынвая и u(0,t) = u(l,t)
u_tx = "x*sin(t)"

alpha = "-1"
beta = "1"
h = float(l) / n
tau = float(T) / m
gamma = tau / h**2

print gamma

def j_area(j_u):
    area = 0
    for i in xrange(n + 1):
        if (i == 0) or (i == n):
            area += j_u[i]*j_u[i]
        elif i % 2 == 0:
            area += 2 * j_u[i]*j_u[i]
        else:
            area += 4 * j_u[i]*j_u[i]
    return h * area / 3

def disc(foo):
    g = np.zeros([m+1, n+1])
    for k in xrange(m+1):
        for i in xrange(n+1):
            g[k, i] = np.array(eval(expr(foo).compile(),
                                      {'x': i*h, 't': k*tau, 'sin': math.sin, 'cos': math.cos, 'exp': math.exp}))
    return g

def direct_sys(g):
    u = np.empty(shape=(m+1, n+1))
    u[0, :] = 0
    u[:, 0] = 0
    u[:, n] = 0
    for k in xrange(0, m):
        for i in xrange(1, n):
            u[k+1, i] = gamma*u[k, i-1] + gamma*u[k, i+1] + tau*g[k, i] - (2*gamma - 1)*u[k, i]
    return u

def con_sys(v):
    psi = np.zeros([m+1, n+1])
    psi[m, :] = v
    psi[:, n] = 0
    psi[:, 0] = 0
    for k in reversed(xrange(1, m+1)):
        for i in xrange(1, n):
            psi[k-1, i] = psi[k, i] - gamma*(2*psi[k, i] - psi[k, i-1] - psi[k, i+1])
    return psi

fff = "sin(x)"

def own_GPM():
    a = disc(alpha)
    b = disc(beta)

    #тест режим
    #f_x = direct_sys(disc(u_tx))[m, :]

    #стандартный режим
    f_x = disc(fff)[m, :]
    u_prev = disc("0")
    u_next = np.zeros([m+1, n+1])
    lmbd = 1.
    st = 0
    while True:
        dj_u = con_sys(direct_sys(u_prev)[m, :] - f_x)
        for k in xrange(m+1):
            for i in xrange(n+1):
                u_next[k, i] = u_prev[k, i] - lmbd * 2 * dj_u[k, i]
                if u_next[k, i] < a[k, i]:
                    u_next[k, i] = a[k, i]
                if u_next[k, i] > b[k, i]:
                    u_next[k, i] = b[k, i]
        if j_area(direct_sys(u_next)[m, :]-f_x) > j_area(direct_sys(u_prev)[m, :]-f_x):
            u_next = copy(u_prev)
            lmbd *= 0.95
            #lmbd *= 1./(st + 1)
        u_prev = copy(u_next)
        st += 1
        print "J: " + str(j_area(direct_sys(u_prev)[m, :]-f_x)), "lambda: " + str(lmbd)
        if st == 1000:
            break
    return u_prev
a = own_GPM()

print direct_sys(disc(u_tx))[m, :]

#тест режим
#print j_area(direct_sys(a)[m, :]-direct_sys(disc(u_tx))[m, :])

#стандартный режим
print j_area(direct_sys(a)[m, :]-disc(fff)[m, :])

print "\n"

pd.DataFrame(data=a, index=[x for x in range(m+1)], columns=[x for x in range(n+1)]).to_csv("U.csv.txt")
pd.DataFrame(data=direct_sys(a), index=[x for x in range(m+1)], columns=[x for x in range(n+1)]).to_csv("direct_U.csv.txt")
pd.DataFrame(data=disc(fff), index=[x for x in range(m+1)], columns=[x for x in range(n+1)]).to_csv("fx.csv.txt")
