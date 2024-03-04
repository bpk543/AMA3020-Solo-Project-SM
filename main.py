from sys import platlibdir
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import math

###################################   Taylor   ######################################


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


z_values = np.linspace(-3, 3, 100)
N_values = [5, 10, 15, 20, 30]

for N in N_values:
    err_values = [errfun(z, N) for z in z_values]
    plt.plot(z_values, err_values, label=f'N={N}')
plt.ylim(-2, 2)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.title('Plot of Taylor series of the Error Function for Different Values of N')
#plt.savefig('Taylor_series_erf2.png', dpi = 500)

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
  
eps=1e-15
eps3= 1e-13
diff01=ErrSequence(0.1,eps3)
diff10=ErrSequence(1.0,eps3)
#diff20=ErrSequence(2.0,eps3)
diff30=ErrSequence(3.0,eps3)
print('this is test:',diff30)
print(len(diff30))
#plt.semilogy(range(1,len(diff01)+1),diff01,'-d', label='z=0.1')
plt.semilogy(range(1,len(diff10)+1),diff10,'-d', label='z=1.0')
#plt.semilogy(range(1,len(diff20)+1),diff20,'-x', label='z=2.0')
plt.semilogy(range(1,len(diff30)+1),diff30,'-d', color = 'r', label='z=3.0')
plt.xlabel('The order of the approximating Taylor Polynomial,N')
plt.ylabel('The absolute difference between $erf(z)$ and $P_N(z)$')
plt.title('The absolute value of the errors in the Taylor Series')
plt.legend()

#plt.savefig('absolute difference.png', dpi=300)

plt.close()


#############################   Asymptotic expansion   ##########################
'''

x = np.linspace(-3, 3, 100)
y = erf(x)

plt.plot(x, y)


def double_factoriall(n):
  mult = 1
  for i in range(1,n+1):
    if i%2 != 0:
      mult = mult*i
  return mult

print(double_factoriall(5))

def error_function_asymptotic(x, terms):
    count = 0
    const = (np.exp(-x**2))/(x*np.sqrt(np.pi))
    for ii in range(terms):
      val = (-1**ii)((double_factoriall(2*ii-1))/(2*(x**2))**ii)
      count += val
    return 1-(const*count)


#result = 1
#for n in range(1, terms + 1):
    #result -= math.exp(-x**2) * ((-1)**n * math.factorial(2*n) / (2**(2*n) * math.factorial(n) * #math.factorial(n) * x**(2*n)))


z_values = np.linspace(-3, 3, 100)
N_values = [5, 10, 15, 20, 30]
#y = math.erf(z_values)
for N in N_values:
    err_values = [error_function_asymptotic(z, N) for z in z_values]
    plt.plot(z_values, err_values, label=f'N={N}')

#plt.plot(z_values, y, label='erf(z)')
plt.xlim(-3,3)
plt.ylim(-2, 2)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.title('Plot of Asymptotic expansion of the Error Function for Different Values of N')
#plt.savefig('Asymptotic expansion_erf.png', dpi = 500)

plt.close()

def ErSequence(z,eps):
  list=[]
  diff = 10.0 
  N = 1       
  while diff > eps:
        diff = abs(math.erf(z)-error_function_asymptotic(z,N))
        list.append(diff)
        N += 1
  return list 
#
# MAIN SCRIPT
#
  
eps=1e-15
eps3= 1e-13
#dif01=ErSequence(0.1,eps3)
#dif10=ErSequence(1.0,eps3)
#diff20=ErrSequence(2.0,eps3)

#dif30=ErSequence(3.0,eps3)
dif10=ErSequence(10.0,eps3)
dif5=ErSequence(5.0,eps3)
#print('this is test:',dif30)
#print(len(diff30))

#plt.semilogy(range(1,len(diff01)+1),diff01,'-d', label='z=0.1')
plt.semilogy(range(1,len(dif10)+1),dif10,'-d', label='z=10.0')
plt.semilogy(range(1,len(dif5)+1),dif5,'-d', label='z=5.0')
#plt.semilogy(range(1,len(diff20)+1),diff20,'-x', label='z=2.0')
#plt.semilogy(range(1,len(dif30)+1),dif30,'-d', color = 'r', label='z=3.0')
plt.xlabel('The order of the approximating Taylor Polynomial,N')
plt.ylabel('The absolute difference between $erf(z)$ and $P_N(z)$')
plt.title('The absolute value of the errors in the Taylor Series')
plt.legend()

plt.savefig('asymtotic_absolute_ difference.png', dpi=300)

plt.close()

'''

################################   Burmann ###########################################



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

plt.close()

'''

def errSequence(z, eps):
    list = []
    diff = 10.0 
    N = 1       
    while diff > eps :
        diff = abs(math.erf(z) - error_function(z, N))
        list.append(diff)
        print(diff)
        N += 1
    return list 

'''
eps1 = 0.3
eps = 1e-15
eps3 = 1e-13
eps4 = 1e-5
#two_diff01 = errSequence(3, eps4)
#two_diff10 = errSequence(1.0, eps4)
#two_diff30 = errSequence(3.0, eps4)

#one_diff30 = errSequence(3.0, eps1)
#print(one_diff30)

#print(two_diff30)
#plt.semilogy(range(1, len(two_diff10) + 1), two_diff10, '-d', label='z=1.0')
#plt.semilogy(range(1, len(two_diff30) + 1), two_diff30, '-d', color='r', label='z=3.0')

plt.title('The absolute value of the errors in the Burmann Series')
#plt.legend()
#plt.savefig('Burmann_absolute_difference.png', dpi=300)
#plt.show()

plt.close()





'''

def error_function_burmann(x):
  const = (2/np.sqrt(np.pi))*np.sqrt(1-np.exp(-x**2))
  val = (np.sqrt(np.pi)) / 2 + (31/200) * np.exp(-x**2) - (341/8000)*np.exp(-2*x**2)
  return (const*val)

print(error_function_burmann(5))
print(erf(5))

x = np.linspace(-3,3,1000)
y = []
for kk in x:
    y.append(abs(erf(x)-error_function_burmann(x)))


plt.plot(x,y)

plt.savefig('test.png')
'''


import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

def error_function_burmann(x):
    const = (2/np.sqrt(np.pi))*np.sqrt(1-np.exp(-x**2))
    val = (np.sqrt(np.pi)) / 2 + (31/200) * np.exp(-x**2) - (341/8000)*np.exp(-2*x**2)
    return (const*val)

x = np.linspace(0, 4, 1000)
y = abs(erf(x) - error_function_burmann(x))

plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('Absolute Difference')
plt.title('Absolute Difference between erf(x) and Practical Burmann(x)')
plt.grid(True)
plt.savefig('Practical Burmann(x).png')

plt.close()

print(max(y))


marks = [1086, 1044, 1221, 3669, 1000, 2843]
# initialize maximum value and index
maximum_val = y[0]
maximum_index = 0
# Initialize loop variable
i = 1
# iterate through the list using while loop
while i < len(y):
    if y[i] > maximum_val:
        maximum_val = y[i]
        maximum_index = i
    i += 1
print(f"Maximum Value: {maximum_val}")
print(f"Maximum Index position: {maximum_index}")
print(x[maximum_index])
print(x[468])


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

x_values = np.linspace(0, 4, 1000)
N_values = [5, 10, 15, 20, 30]

for N in N_values:
    err_values = [abs(erf(x) - errfun(x, N)) for x in x_values]
    plt.plot(x_values, err_values, label=f'N={N}')

plt.ylim(0, 0.13)
plt.xlabel('x')
plt.ylabel('Absolute Difference')
plt.title('Absolute Difference between erf(x) and the Taylor series for Different N values')
plt.legend()
plt.grid(True)
plt.savefig('func.png')

plt.close()

x_value = 1.4054054
N_values = [5, 10, 15, 20, 30]

print(f"x = {x_value}")
for N in N_values:
    y_value = abs(erf(x_value) - errfun(x_value, N))
    print(f"Taylor N={N}, y={y_value}")

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.special import erf

def errfun(z, N):
    summation = 0
    for n in range(N):
        a = (-1) ** n
        b = z ** ((2 * n) + 1)
        c = math.factorial(n)
        d = (2 * n) + 1
        P = (a * b) / (c * d)
        summation += (2 / np.sqrt(np.pi)) * P
    return summation

x_values = np.linspace(0, 4, 1000)
N_values = [5, 10, 15, 20, 30]

for N in N_values:
    err_values = [abs(erf(x) - errfun(x, N)) for x in x_values]
    plt.plot(x_values, err_values, label=f'N={N}')

plt.ylim(0, 40)  # Set y-axis limits from None to 40
plt.xlabel('x')
plt.ylabel('Absolute Difference')
plt.title('Absolute Difference between erf(x) and Taylor Series for Different N values')
plt.legend()
plt.grid(True)
plt.savefig('erf(x) and Taylor Series.png')
plt.close()


'''
def error_function_asymptotic(x, terms):
  result = 1
  for n in range(1, terms + 1):
      result -= math.exp(-x**2) * ((-1)**n * math.factorial(2*n) / (2**(2*n) * math.factorial(n) * math.factorial(n) * x**(2*n)))
  return result

x_values = np.linspace(0, 10, 1000)
N_values = [5, 10, 15, 20, 30]

for N in N_values:
    err_values = [abs(erf(x) - error_function_asymptotic(x, N)) for x in x_values]
    plt.plot(x_values, err_values, label=f'N={N}')

plt.xlabel('x')
plt.ylabel('Absolute Difference')
plt.title('Absolute Difference between erf(x) and error_function_asymptotic for Different N values')
plt.legend()
plt.grid(True)
plt.savefig('test3.png')

plt.close()
'''


'''
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.special import erf

def double_factoriall(n):
  mult = 1
  for i in range(1,n+1):
    if i%2 != 0:
      mult = mult*i
  return mult

def error_function_asymptotic_sam(x, terms):
  count = 0
  const = (np.exp(-x**2))/(x*np.sqrt(np.pi))
  for ii in range(terms):
    val = (-1**ii)((double_factoriall(2*ii-1))/(2*(x**2))**ii)
    count += val
  return 1-(const*count)

def error_function_asymptotic(x, terms):
    result = 1
    for n in range(1, terms + 1):
        denominator = 2**(2*n) * math.factorial(n) * math.factorial(n) * x**(2*n)
        if denominator != 0:
            result -= math.exp(-x**2) * ((-1)**n * math.factorial(2*n) / denominator)
    return result

x_values = np.linspace(0, 0.1, 1000)
N_values = [5, 10, 15, 20, 30]

for N in N_values:
    err_values = [abs(erf(x) - error_function_asymptotic(x, N)) for x in x_values]
    plt.plot(x_values, err_values, label=f'N={N}')

plt.xlabel('x')
plt.ylabel('Absolute Difference')
plt.title('Absolute Difference between erf(x) and error_function_asymptotic for Different N values')
plt.legend()
plt.grid(True)
plt.savefig('test3.png')
plt.close()

#########################################   Sam   ###################################

def double_factoriall(n):
  mult = 1
  for i in range(1,n+1):
    if i%2 != 0:
      mult = mult*i
  return mult

def error_function_asymptotic_sam(x, terms):
  count = 0
  const = (np.exp(-x**2))/(x*np.sqrt(np.pi))
  for ii in range(terms):
    val = (-1**ii)*((double_factoriall(2*ii-1))/(2*(x**2))**ii)
    count += val
  return 1-(const*count)

def error_function_asymptotic_sam2(x, terms):
  count = 0
  const = (np.exp(-x**2))/(x*np.sqrt(np.pi))
  for ii in range(terms):
      denominator = 2 * (x ** 2)
      if denominator != 0:
          val = ((-1) ** ii) * ((double_factoriall(2*ii-1)) / denominator ** ii)
          count += val
  return 1 - (const * count)

x_values = np.linspace(0, 0.1, 1000)
N_values = [5, 10, 15, 20, 30]

for N in N_values:
    err_values = [abs(erf(x) - error_function_asymptotic_sam2(x, N)) for x in x_values]
    plt.plot(x_values, err_values, label=f'N={N}')

plt.xlabel('x')
plt.ylabel('Absolute Difference')
plt.title('Absolute Difference between erf(x) and error_function_asymptotic for Different N values')
plt.legend()
plt.grid(True)
plt.savefig('test_sam2.png')
plt.close()

'''

##############################    Burmann absolte   #################


def error_function(z, N):
  result = 0
  for n in range(N):
      coefficient = (-1)**n / (math.factorial(n) * (2 * n + 1))
      term = coefficient * (z / math.sqrt(2))**(2 * n + 1)
      result += term
  return 2 / math.sqrt(math.pi) * result

x_values = np.linspace(0, 4, 1000)
N_values = [5, 10, 15, 20, 30]

for N in N_values:
    err_values = [abs(erf(x) - error_function(x, N)) for x in x_values]
    plt.plot(x_values, err_values, label=f'N={N}')

plt.ylim(None, 40)
plt.xlabel('x')
plt.ylabel('Absolute Difference')
plt.title('Absolute Difference between erf(x) and Burmann Series for Different N values')
plt.legend()
plt.grid(True)
plt.savefig('erf(x) and Burmann Series.png')

plt.close()

x_value = 1.4054054
N_values = [4, 5, 10, 15, 20, 30]

print(error_function(x_value, 4))
print(error_function(x_value, 5))
print(erf(x_value))
print(error_function_burmann(x_value))

print(f"x = {x_value}")
for N in N_values:
    y_value = abs(erf(x_value) - error_function(x_value, N))
    print(f"Burmann N={N}, y={y_value}")


print(f"x = {x_value}")
for N in N_values:
    y_value = errfun(x_value, N)
    print(f"Taylor 1.4054 N={N}, y={y_value}")

print(f"x = {x_value}")
for N in N_values:
    y_value = error_function(x_value, N)
    print(f"Burmann 1.4054 N={N}, y={y_value}")

print('real erf(1.4054) = ', erf(1.4054))