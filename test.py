from scipy.optimize import minimize

def fun(x):
    return 1*x**2 + 2*x + 3

x = minimize(fun, 10)

print(x)

##