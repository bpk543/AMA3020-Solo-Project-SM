import numpy as np
from sys import platlibdir

import matplotlib.pyplot as plt
from scipy.special import erf
import math

############################    Burmann   ################################


def burmann1(x):
  const = (2/np.sqrt(np.pi))*np.sqrt(1-np.exp(-x**2))
  val = 1
  return (const*val)

def burmann2(x):
  const = (2/np.sqrt(np.pi))*np.sqrt(1-np.exp(-x**2))
  val = 1-(1/12)*(1-np.exp(-x**2))
  return (const*val)

def burmann3(x):
  const = (2/np.sqrt(np.pi))*np.sqrt(1-np.exp(-x**2))
  val = 1-1/12*(1-np.exp(-x**2)) -(7/480)*(1-np.exp(-x**2))**2
  return (const*val)

def burmann4(x):
  const = (2/np.sqrt(np.pi))*np.sqrt(1-np.exp(-x**2))
  val = 1-1/12*(1-np.exp(-x**2)) -(7/480)*(1-np.exp(-x**2))**2 -(5/896)*(1-np.exp(-x**2))**3
  return (const*val)

def burmann5(x):
  const = (2/np.sqrt(np.pi))*np.sqrt(1-np.exp(-x**2))
  val = 1-1/12*(1-np.exp(-x**2)) -(7/480)*(1-np.exp(-x**2))**2 -(5/896)*(1-np.exp(-x**2))**3 -(787/276840)*(1-np.exp(-x**2))**4
  return (const*val)

print(burmann5(0.5))
print(burmann5(-0.5))
print(-burmann5(-0.5))



x = np.linspace(0, 3, 1000)
y = np.linspace(-3,0,1000)
plt.plot(x, burmann1(x), color='#1f77b4', label="N=1")
plt.plot(x, burmann2(x), color='#ff7f0e', label="N=2")
plt.plot(x, burmann3(x), color='#2ca02c', label="N=3")
plt.plot(x, burmann4(x), color='#d62728', label="N=4")
plt.plot(x, burmann5(x), color='#9467bd', label="N=5")
plt.plot(y, -burmann1(y), color='#1f77b4')
plt.plot(y, -burmann2(y), color='#ff7f0e')
plt.plot(y, -burmann3(y), color='#2ca02c')
plt.plot(y, -burmann4(y), color='#d62728')
plt.plot(y, -burmann5(y), color='#9467bd')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Plot of Burmann Series of the Error Function for Different Values of N')
plt.savefig('burmann.png')

plt.close()

  

x = np.linspace(0, 4, 1000)
a = abs(erf(x) - burmann1(x))
b = abs(erf(x) - burmann2(x))
c = abs(erf(x) - burmann3(x))
d = abs(erf(x) - burmann4(x))
e = abs(erf(x) - burmann5(x))

print(e)


#plt.plot(x, y)
plt.plot(x, a, color='#1f77b4', label="N=1")
plt.plot(x, b, color='#ff7f0e', label="N=2")
plt.plot(x, c, color='#2ca02c', label="N=3")
plt.plot(x, d, color='#d62728', label="N=4")
plt.plot(x, e, color='#9467bd', label="N=5")
plt.legend()
plt.xlabel('x')
plt.ylabel('Absolute Difference')
plt.title('Absolute Difference between erf(x) and Burmann Series for Different N values')
plt.grid(True)
plt.savefig('absolute_Burmann(x)2.png')

plt.close()

plt.plot(x, a, color='#1f77b4', label="N=1")
plt.plot(x, b, color='#ff7f0e', label="N=2")
plt.plot(x, c, color='#2ca02c', label="N=3")
plt.plot(x, d, color='#d62728', label="N=4")
plt.plot(x, e, color='#9467bd', label="N=5")
plt.legend()
plt.ylim(None, 40)
plt.xlabel('x')
plt.ylabel('Absolute Difference')
plt.title('Absolute Difference between erf(x) and Burmann Series')
plt.grid(True)
plt.savefig('absolute_Burmann(x)40.png')

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

x = np.linspace(0, 5, 1000)
y = np.linspace(-1,5,1000)
q = np.linspace(-1,0,1000)
plt.plot(x, burmann5(x), 'r--', label="Burmann (N=5)")
plt.plot(q, -burmann5(q), 'r--')
plt.plot(y, errfun(y,5), 'b--', label='Taylor (N=5)')
plt.plot(y, erf(y), 'k', label='erf(x)')
plt.legend()
plt.ylim(-1.1, 2)
plt.xlim(-1,None)
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Expansion of erf(x) at x=0')
#plt.grid(True)
plt.savefig('comparison4.png')

plt.close()

def practical(x):
  return (2/np.sqrt(np.pi))*np.sqrt(1-np.exp(-x**2))*((np.sqrt(np.pi))/2 + (31/200)*np.exp(-x**2) - (341/8000)*np.exp(-2*x**2) )

x = np.linspace(0, 5, 1000)
plt.plot(x, abs(erf(x) - practical(x)))
plt.xlabel('x')
plt.ylabel('y')
plt.title('Absolute Difference between erf(x) and practical Burmann')
#plt.grid(True)
plt.savefig('practical.png')

