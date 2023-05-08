from sympy import *
from sympy.abc import x, y

def quasi_polydivision(f,g):
    # generates k, h with deg(LC(f)^d * g + h*f) < deg(f) 
    if f == 0:
        raise ValueError('No devision by 0 possible.') 
    m = degree(f, gen=x)
    n = degree(g, gen=x)
    if n < m:
        return 0,0
    c = LC(g)
    d = LC(f)
    l = n-m
    k,h = quasi_polydivision(f,d*g-c*f*x**l)
    return k+1, h-d**k*c*x**l

def new_prime(d):
    # Returns the smallest prime number which is no divisor of d.
    if d == 0:
        raise ValueError('No devisions by 0 possible.')
    i = 1
    while True:
        p = prime(i)
        if d % p == 0:
            i += 1
            continue
        else:
            return p

def coefficent_polynomial(m,i):
    # generates g_i such that (1+m)^i = 1+m*g_i
    if i == 0:
        return 0
    else:
        g = coefficent_polynomial(m,i-1)
        return m*g+g+1

# In the following M,nu is an input and the output has the form  "True, *x" or "False, *x".
# If the output has the form "True, *x", the existence statement is proved and *x is a widness of the existence.
# If the output has the form "False, *x" is means that (M,nu) is not an explicit maximal ideal and *x is a widness for it.
# Here we have the following cases:
    # False, 0, 0      =>  0 not in M
    # False, 0, 1      =>  1 in M
    # False, 1, a, b   =>  a,b in M but a+b not in M
    # False, 2, a, b   =>  a in M but ab not in M
    # False, 3, a      =>  a not in M and a*nu(a)-1 not in M
    
def maxideal_primeideal(factors, M, nu):
    # takes a list factors of elements in Z[X] and if its product is in M, it either return an element in factors which is also in M, 
    # or it return evidence that M,nu is not an explicit maximal ideal.
    if factors is []:
        return False, 0,1
    a = factors[0]
    b = 1
    for i in factors[1:]:
        b *= i
    if M(a):
        return True, a
    if M(b):
        return maxideal_primeideal(factors[1:], M, nu)
    if not M(a*nu(a)*b):
        return False, 2, ab, nu(a)    
    if not M(a*nu(a)-1):
        return False, 3,a
    if not M(b-a*nu(a)*b):
        return False, 2, a*nu(a)-1, -b
    return False, 1, a*nu(a)*b, b-a*nu(a)*b

def MaxZX(M, nu):
    # takes M, nu and either return a prime number p in M or evidence that (M,nu) is not an explicit maximal ideal. 
    
    if M(x):
        f = x
    else:
        f = x*nu(x) - 1
        if not M(f):
            return False, 3, x
    d = LC(Poly(f,x))
    q = new_prime(d)
    if M(q):
        return True, q
    m = q*nu(q)-1
    if not M(m):
        return False, 3, q
    n = degree(f)
    K = []
    H = []
    A = []
    for i in range(n):
        k_i, h_i = quasi_polydivision(f,nu(q)*x**i)
        K.append(k_i)
        H.append(-h_i)
        F = d**k_i*nu(q)*x**i+h_i*f
        coeffs = Poly(F,x).all_coeffs()
        coeffs.reverse()
        line_i = coeffs + [0] * (n-len(coeffs))
        line_i[i] -= d**k_i*y
        A.append([-e for e in line_i])
    B = Matrix(A).adjugate()
    factor_f = (B*Matrix(H)).subs(y, nu(q))[0]*q**n
    char_poly = Matrix(A).det()
    factor_m = 0
    Bs = Poly(char_poly, y).all_coeffs()
    Bs.reverse()
    for i in range(n+1):
        factor_m -= q**(n-i)*Bs[i]*coefficent_polynomial(m,i)
    if  not M(factor_f * f ):
        return False, 2, f, factor_f
    if not M(factor_m *m):
        return False, 2, m, factor_m
    if not M(factor_f * f +  factor_m *m):
        return False, 1, factor_f * f,  factor_m *m
    if factor_f * f +  factor_m *m == 1 or factor_f * f +  factor_m *m == -1:
        return False, 0
    return maxideal_primeideal(factorint(simplify(factor_f * f +  factor_m *m),multiple = True), M, nu)

############################## EXAMPLES ##############################

def M(f):
    k,h = quasi_polydivision(Poly(x**2+1,x),f)
    g = f+h*(x**2+1)
    L = poly(g, x).all_coeffs()
    if len(L) == 0:
        return True
    if len(L) == 1:
        if L[0] %3 == 0:
            return True
        else:
            return False
    if L[0] % 3 == 0 and L[1] % 3 == 0:
        return True
    else:
        return False
    
def nu(f):
    k,h = quasi_polydivision(Poly(x**2+1,x),f)
    g = f+h*(x**2+1)
    L = poly(g, x).all_coeffs()
    L = [a % 3 for a in L]
    if len(L) == 0:
        return 0
    if len(L) == 1:
        if L[0] == 1:
         return x**2 + 2
        if L[0] == 2:
            return x**2 
    if L[0] == 1:
        if L[1] == 0:
            return - x
        if L[1] == 1:
            return x - 1
        else:
            return x + 1
    else:
        if L[1] == 0:
            return x
        if L[1] == 1:
            return -x - 1
        else:
            return -x + 1
    
print(MaxZX(M, nu))
