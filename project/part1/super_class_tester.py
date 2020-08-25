from super_class import integral
from numpy import sqrt, pi, exp

sigma1 = 1
sigma2 = 1
sigma3 = 1

a1 = -1
a2 = -2
a3 = -3

b1 = -a1
b2 = -a2
br = -a3

mu = 0

test1 = integral(sigma1, mu, a1, b1)
int_test1 = test1.P()

print(int_test1)

test2 = integral(sigma2, mu, -sigma2, sigma2)
int_test2 = test2.P()

print(int_test2)

test3 = integral(sigma2, mu, -sigma3, sigma3)
int_test3 = test3.P()

print(int_test3)
