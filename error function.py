import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import math 

x = np.linspace(-3, 3, 100)
y = erf(x)

plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('Error Function')
plt.title('Error Function Plot')
#plt.savefig('erf_plot.png')
plt.close()

x = np.linspace(-3, 3, 100)
y = erf(x)
z = 2/np.pi * np.arctan(2*x*(1+x**4))
plt.plot(x, y, label='erf(x)')
#plt.plot(x, z, label='2/pi * arctan(2*x*(1+x**4))')
plt.plot(x, z, label=r'$\frac{2}{\pi}  \arctan(2x(1+x^4))$')
plt.xlabel('x')
plt.ylabel('y')
plt.legend(fontsize="13")
plt.title('Error Function vs Approximation Plot')
#plt.savefig('erf+approx2.png')
plt.close()

def errfun(z,N):
  sum=0
  for n in range (0,N):
    a=(-1)**n
    b=z**((2*n)+1)
    c=math.factorial(n)
    d=(2*n)+1
    P=(a*b)/(c*d)
    sum+=(2/math.sqrt(math.pi))*P

  return sum


print( errfun(1,9) , "= approximate")
print(erf(1) , "= real")


x = np.linspace(-3, 3, 100)
y = erf(x)

plt.plot(x, y)


z_values = np.linspace(-3, 3, 100)
N_values = [14, 42]

for N in N_values:
    err_values = [errfun(z, N) for z in z_values]
    plt.plot(z_values, err_values, label=f'N={N}')
plt.ylim(-2, 2)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.title('Plot of Taylor series of the Error Function for Different Values of N')
plt.savefig('Taylor_series_erf2.png', dpi = 500)

'''
z_values = np.linspace(-3, 3, 100)
N_values = [5, 10, 15, 20]

for N in N_values:
    err_values = [errfun(z, N) for z in z_values]
    plt.plot(z_values, err_values, label=f'N={N}')
plt.ylim(-2, 2)
plt.xlabel('z')
plt.ylabel('errfun(z, N)')
plt.legend()
plt.title('Plot of errfun for Different Values of N')

plt.savefig('errfun_plot.png')
'''
plt.close()

#  
def ErrSequence(z,eps):
  list=[]
  diff = 10.0 
  N = 1       
  while diff > eps:
        diff = abs(math.erf(z)-errfun(z,N))
        list.append(diff)
        N += 1
  return list 
#
# MAIN SCRIPT
#
import matplotlib.pyplot as pl  
eps=1e-15
eps3= 1e-13
diff01=ErrSequence(0.1,eps3)
diff10=ErrSequence(1.0,eps3)
#diff20=ErrSequence(2.0,eps3)
diff30=ErrSequence(3.0,eps3)
print('this is test:',diff30)
print(len(diff30))
#pl.semilogy(range(1,len(diff01)+1),diff01,'-d', label='z=0.1')
pl.semilogy(range(1,len(diff10)+1),diff10,'-d', label='z=1.0')
#pl.semilogy(range(1,len(diff20)+1),diff20,'-x', label='z=2.0')
pl.semilogy(range(1,len(diff30)+1),diff30,'-d', color = 'r', label='z=3.0')
pl.xlabel('The order of the approximating Taylor Polynomial,N')
pl.ylabel('The absolute difference between $erf(z)$ and $P_N(z)$')
pl.title('The absolute value of the errors in the Taylor Series')
pl.legend()

#pl.savefig('absolute difference.png', dpi=300)


#print(ErrSequence(0.1,1e-13)[5])
print(ErrSequence(1.0,1e-13)[12])
print(ErrSequence(1.0,1e-13)[13])
print(ErrSequence(1.0,1e-13)[14])

#print(ErrSequence(2.0,1e-15)[5])
#print(ErrSequence(3.0,1e-13)[5])
pl.close()
#print(ErrSequence(1.0,1e-15)[4])
#print(ErrSequence(1.0,1e-15)[14])
print(ErrSequence(3.0,1e-13)[40])
print(ErrSequence(3.0,1e-13)[41])
print(ErrSequence(3.0,1e-13)[42])

#print(ErrSequence(2.0,1e-15)[28])

#########################################################################################
################  gauss_hermite_quadrature  #############################################
#########################################################################################
'''

import numpy as np
from scipy.special import erf

def test(x):
  return 1/((x-1/2)**2+4)

def gauss_hermite_quadrature(func, num_points):
    # Calculate the nodes and weights for Gauss-Hermite quadrature
    nodes, weights = np.polynomial.hermite.hermgauss(num_points)

    # Perform the numerical integration using Gauss-Hermite quadrature
    result = np.sum(weights * func(nodes))

    return result

# Define the error function
def error_function(x):
    return erf(x)

def one(x):
  return 1

# Set the number of points for Gauss-Hermite quadrature
num_points = 10

# Perform the numerical integration of the error function

################################ test
test_result2 = gauss_hermite_quadrature(test, 2)
test_result3 = gauss_hermite_quadrature(test, 3)

print(f'Test Result for num_points=2: {test_result2}')
print(f'Test Result for num_points=3: {test_result3}')

############################# one 

func_1_result2 = gauss_hermite_quadrature(one, 2)
func_1_result3 = gauss_hermite_quadrature(one, 3)
func_1_result10 = gauss_hermite_quadrature(one, 10)

print(f'func_1 Result for num_points=2: {func_1_result2}')
print(f'func_1 Result for num_points=3: {func_1_result3}')
print(f'func_1 Result for num_points=10: {func_1_result10}')


##############################################

integral_result = gauss_hermite_quadrature(error_function, 10)
integral_result3 = gauss_hermite_quadrature(error_function, 3)
print(f'Numerical Integration Result for num_points=3: {integral_result3}')
integral_result2 = gauss_hermite_quadrature(error_function, 2)
print(f'Numerical Integration Result for num_points=2: {integral_result2}')

integral_result9 = gauss_hermite_quadrature(error_function, 9)
print(f'Numerical Integration Result for num_points=9: {integral_result9}')


integral_result100 = gauss_hermite_quadrature(error_function, 100)
print(f'Numerical Integration Result for num_points=100: {integral_result100}')

print(f'Numerical integration result using Gauss-Hermite quadrature: {integral_result}')


# Plotting the Gauss-Hermite quadrature for different num_points
num_points_list = [5, 10, 15, 20, 25, 30]
integral_results = []

for num_points in num_points_list:
    integral_result = gauss_hermite_quadrature(error_function, num_points)
    integral_results.append(integral_result)

plt.plot(num_points_list, integral_results, marker='o', linestyle='-', color='b')
plt.xlabel('Number of Points')
plt.ylabel('Integral Result')
plt.title('Gauss-Hermite Quadrature for Different Number of Points')
plt.grid(True)
plt.savefig('gauss_hermite_quadrature.png')

################################################################################

import numpy as np
from scipy.special import erf
from scipy.integrate import quad

# Define the integrand function
def integrand(t):
    return np.exp(-t**2)

# Number of quadrature points
n = 100

# Compute the nodes (x) and weights (w) for Gauss-Hermite quadrature
x, w = np.polynomial.hermite.hermgauss(n)

# Scale the nodes and weights to the interval [0, inf)
x_scaled = np.sqrt(2) * x
w_scaled = w / np.sqrt(np.pi)

# Evaluate the integral using Gauss-Hermite quadrature
result_approx = np.dot(w_scaled, integrand(x_scaled))

# Evaluate the error function using scipy's erf function
result_exact = erf(x_scaled[-1])

print("Approximate result using Gauss-Hermite quadrature:", result_approx)
print("Exact result using scipy's erf function:", result_exact)

# Plotting the Gauss-Hermite quadrature for different number of quadrature points vs the actual error function
num_points_list = [5, 10, 15, 20, 25, 30]
integral_results = []

for num_points in num_points_list:
    integral_result = gauss_hermite_quadrature(error_function, num_points)
    integral_results.append(integral_result)

plt.plot(num_points_list, integral_results, marker='o', linestyle='-', color='b')
plt.xlabel('Number of Quadrature Points')
plt.ylabel('Integral Result')
plt.title('Gauss-Hermite Quadrature vs Actual Error Function')
plt.grid(True)
plt.savefig('gauss_hermite_quadrature_vs_erf.png')

plt.close()
# Plotting the Gauss-Hermite quadrature for different number of quadrature points vs the actual error function
num_points_list = [100]
x_values = np.linspace(-3, 3, 100)

for num_points in num_points_list:
    integral_results = []
    for x in x_values:
        integral_result = gauss_hermite_quadrature(lambda t: np.exp(-t**2) * erf(x*t), num_points)
        integral_results.append(integral_result)
    plt.plot(x_values, integral_results, label=f'Number of Quadrature Points = {num_points}')

plt.xlabel('x')
plt.ylabel('Integral Result')
plt.title('Gauss-Hermite Quadrature for Different Number of Quadrature Points')
plt.grid(True)
plt.legend()
plt.savefig('gauss_hermite_quadrature_vs_erf_x_range.png')

plt.close()

###############################################################################


import numpy as np

def ff(x):
    return np.cos(x)

def one(x):
  return 1

def g(x):
    return np.exp(-x**2) * ff(x)

# Interval for Chebfun
interval_chebfun = [-np.inf, np.inf]

# Interval for limited integration
interval_limited = [-6, 6]

# Real interval
interval_real = [0, 3]

# Exact value
exact = np.sqrt(np.pi) * np.exp(-1/4)

# Chebfun integration
integral_chebfun = np.sum(g(np.linspace(interval_chebfun[0], interval_chebfun[1], 1000)))

# Limited interval integration
integral_limited = np.sum(g(np.linspace(interval_limited[0], interval_limited[1], 1000)))

integral_real = np.sum(g(np.linspace(interval_real[0], interval_real[1], 1000)))

print("Chebfun Integration (infinite interval):", integral_chebfun)
print("Chebfun Integration (limited interval):", integral_limited)
print("Chebfun Integration (real interval):", integral_real)

print("\nHermite Quadrature Errors:")
print("    n        error")
for n in range(1, 13):
    s, w = np.polynomial.hermite.hermgauss(n)
    In = np.sum(w * ff(s))
    err = In - exact
    print(f'{n:3d} {err:19.15f}')

'''
########################################  Burmann  #############################
import math

def error_function(z, N):
    result = 0
    for n in range(N):
        coefficient = (-1)**n / (math.factorial(n) * (2 * n + 1))
        term = coefficient * (x / math.sqrt(2))**(2 * n + 1)
        result += term
    return 2 / math.sqrt(math.pi) * result

# Example usage:
x = 1.5
terms = 50
print("erf(", x, ") =", error_function(z,N))


z_values = np.linspace(-3, 3, 100)
N_values = [5, 10, 15, 20, 30]

for N in N_values:
    err_values = [error_function(z, N) for z in z_values]
    plt.plot(z_values, err_values, label=f'N={N}')
plt.ylim(-2, 2)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.title('Plot of Burmann series of the Error Function for Different Values of N')
#plt.savefig('Burmann_series_erf.png', dpi = 500)

plt.close()

'''
def errSequence(z,eps):
  list=[]
  diff = 10.0 
  N = 1       
  while diff > eps:
        diff = abs(math.erf(z)-error_function(z,N))
        list.append(diff)
        N += 1
  return list 
#
# MAIN SCRIPT
#
import matplotlib.pyplot as pl  
eps=1e-15
eps3= 1e-13
eps4 = 1e-5
diff01=errSequence(0.1,eps4)
diff10=errSequence(1.0,eps4)
#diff20=ErrSequence(2.0,eps3)
diff30=errSequence(3.0,eps4)
print(diff30)
#pl.semilogy(range(1,len(diff01)+1),diff01,'-d', label='z=0.1')
pl.semilogy(range(1,len(diff10)+1),diff10,'-d', label='z=1.0')
#pl.semilogy(range(1,len(diff20)+1),diff20,'-x', label='z=2.0')
pl.semilogy(range(1,len(diff30)+1),diff30,'-d', color = 'r', label='z=3.0')
#pl.xlabel('The order of the approximating Taylor Polynomial,N')
#pl.ylabel('The absolute difference between $erf(z)$ and $P_N(z)$')
pl.title('The absolute value of the errors in the Burmann Series')
pl.legend()
pl.savefig('Burmann absolute difference .png', dpi=300)

plt.close()

'''

import math
import numpy as np
import matplotlib.pyplot as plt

def error_function(z, N):
    result = 0
    for n in range(N):
        coefficient = (-1)**n / (math.factorial(n) * (2 * n + 1))
        term = coefficient * (z / math.sqrt(2))**(2 * n + 1)
        result += term
    return 2 / math.sqrt(math.pi) * result

# Example usage:
x = 1.5
N = 50
print("erf(", x, ") =", error_function(x, N))


z_values = np.linspace(-3, 3, 100)
N_values = [5, 10, 15, 20, 30]

for N in N_values:
    err_values = [error_function(z, N) for z in z_values]
    plt.plot(z_values, err_values, label=f'N={N}')
plt.ylim(-2, 2)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.title('Plot of Burmann series of the Error Function for Different Values of N')
#plt.savefig('Burmann_series_erf.png', dpi=500)
#plt.show()

def errSequence(z, eps):
    list = []
    diff = 10.0 
    N = 1       
    while diff > eps:
        diff = abs(math.erf(z) - error_function(z, N))
        list.append(diff)
        N += 1
    return list 

eps = 1e-15
eps3 = 1e-13
eps4 = 1e-5
two_diff01 = errSequence(0.1, eps4)
two_diff10 = errSequence(1.0, eps4)
two_diff30 = errSequence(3.0, eps4)
print(two_diff30)
plt.semilogy(range(1, len(two_diff10) + 1), two_diff10, '-d', label='z=1.0')
plt.semilogy(range(1, len(two_diff30) + 1), two_diff30, '-d', color='r', label='z=3.0')

plt.title('The absolute value of the errors in the Burmann Series')
plt.legend()
#plt.savefig('Burmann_absolute_difference.png', dpi=300)
plt.show()


############################   Asymptotic expansion   #############################
'''
import math

def error_function_asymptotic(x, terms):
    result = 1
    for n in range(1, terms + 1):
        result -= math.exp(-x**2) * ((-1)**n * math.factorial(2*n) / (2**(2*n) * math.factorial(n) * math.factorial(n) * x**(2*n)))
    return result

# Example usage:
x = 10.0
terms = 4  # Increase the number of terms for more accuracy
print("erf(", x, ") =", error_function_asymptotic(x, terms))


z_values = np.linspace(-3, 3, 100)
N_values = [5, 10, 15, 20, 30]
y = math.erf(z_values)
for N in N_values:
    err_values = [error_function_asymptotic(z, N) for z in z_values]
    #plt.plot(z_values, err_values, label=f'N={N}')

plt.plot(z_values, y, label='erf(z)')
plt.xlim(-3,3)
plt.ylim(-2, 2)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.title('Plot of Asymptotic expansion of the Error Function for Different Values of N')
plt.savefig('Asymptotic expansion_erf.png', dpi = 500)

'''

import math
import numpy as np
import matplotlib.pyplot as plt

def error_function_asymptotic(x, terms=3):
    result = 1
    for n in range(1, terms + 1):
        result -= math.exp(-x**2) * ((-1)**n * math.factorial(2*n) / (2**(2*n) * math.factorial(n) * math.factorial(n) * x**(2*n)))
    return result

z_values = np.linspace(-3, 3, 100)
N_values = [5, 10, 15, 20, 30]
y = math.erf(z_values)
for N in N_values:
    err_values = [error_function_asymptotic(z, N) for z in z_values]
    plt.plot(z_values, err_values, label=f'N={N}')

plt.plot(z_values, y, label='erf(z)')
plt.xlim(-3, 3)
plt.ylim(-2, 2)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.title('Plot of Asymptotic Expansion of the Error Function for Different Values of N')
plt.savefig('Asymptotic_expansion_erf2.png', dpi=500)
#plt.show()
