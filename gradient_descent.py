######################################################################
# Metodo de Steepest Descent
# Bruno Nunes Cabral Tenorio
# Metodo para achar um minimo para uma funcao g(x,y,...)
# Maio de 2016
######################################################################
## VOCE VAI PRECISAR DA BIBLIOTECA "sympy"
# http://docs.sympy.org/dev/install.html ou pelo synaptic
from sympy import *
x1,x2,x3,x4,x5 = symbols('x1 x2 x3 x4 x5')
from mpmath import *
#from sympy.mpmath import *
mp.dps = 50; mp.pretty = True
#####################################################################

# definindo a funcao g(x) e seu gradiente
# Se for modificar a funcao e suas variaveis vai ter de modificar em:
# g(x), ng([a,b,...]), diff(lambda x1,...: g(x))

# var = numero de varieveis do sistema
var=raw_input('numero de variaveis do sistema: ' )
var=int(var)
#epsilon = precisao
epsilon = 10**-30
#########################################################
# funcao g(x)
#g(x)=(x1**3 +(x1**2)*x2 -x1*x3+6)**2 +(exp(x1)+exp(x2)-x3)**2 +(x2**2 -2*x1*x3 -4)**2 
str_exp= raw_input('funcao nas variaveis x1,..,xn: ' ) 
expression=sympify(str_exp)
if var ==2:
    def g(X,Y):
        return expression.subs([(x1,X),(x2,Y)])
elif var ==1:
    def g(X):
        return expression.subs([(x1,X)])
elif var ==3:
    def g(X,Y,Z):
        return expression.subs([(x1,X),(x2,Y),(x3,Z)])
elif var == 4:
    def g(X,Y,Z,D):
        return expression.subs([(x1,X),(x2,Y),(x3,Z),(x4,D)])
elif var == 5:
    def g(X,Y,Z,D,E):
        return expression.subs([(x1,X),(x2,Y),(x3,Z),(x4,D),(x5,E)])
# define a funcao que calcula no ponto
def ng(ponto):
    if var == 2:
        return g(x1,x2).subs([(x1,ponto[0]),(x2,ponto[1])])
    elif var ==1:
        return g(x1).subs([(x1,ponto[0])])
    elif var ==3:
        return g(x1,x2,x3).subs([(x1,ponto[0]),(x2,ponto[1]),(x3,ponto[2])])
    elif var ==4:
        return g(x1,x2,x3,x4).subs([(x1,ponto[0]),(x2,ponto[1]),(x3,ponto[2]),(x4,ponto[3])])
    elif var==5:
        return g(x1,x2,x3,x4,x5).subs([(x1,ponto[0]),(x2,ponto[1]),(x3,ponto[2]),(x4,ponto[3]),(x5,ponto[4])])
#########################################################
# calcula a funcao numerica g(lista)
if var == 2:
# (modulo)gradiente da funcao g(x) no ponto (x1,x2)
    gradii=[[0 for m in range(var)] for n in range(var)]
    for j in range(0,var):
        gradii[j][j]=1
    def grad(ponto):
        grad=list()
        modular=list()
        tuple(ponto)
        for elemento in gradii:
            grad.append(diff(lambda x1,x2: g(x1,x2),ponto,tuple(elemento)))
        for element in grad:
            modular.append(element**2)
        modulo =  sum(modular)
        for i in range(0,var):
            grad[i]=grad[i]/sqrt( modulo)
        return grad
elif var == 1:
# (modulo)gradiente da funcao g(x) no ponto (x1)
    gradii=[[0 for m in range(var)] for n in range(var)]
    for j in range(0,var):
        gradii[j][j]=1
    def grad(ponto):
        grad=list()
        modular=list()
        tuple(ponto)
        for elemento in gradii:
            grad.append(diff(lambda x1: g(x1),ponto,tuple(elemento)))
        for element in grad:
            modular.append(element**2)
        modulo =  sum(modular)
        for i in range(0,var):
            grad[i]=grad[i]/sqrt( modulo)
        return grad
elif var == 3:
# (modulo)gradiente da funcao g(x) no ponto (x1,x2,x3)
    gradii=[[0 for m in range(var)] for n in range(var)]
    for j in range(0,var):
        gradii[j][j]=1
    def grad(ponto):
        grad=list()
        modular=list()
        tuple(ponto)
        for elemento in gradii:
            grad.append(diff(lambda x1,x2,x3: g(x1,x2,x3),ponto,tuple(elemento)))
        for element in grad:
            modular.append(element**2)
        modulo =  sum(modular)
        for i in range(0,var):
            grad[i]=grad[i]/sqrt( modulo)
        return grad
elif var == 4:
    gradii=[[0 for m in range(var)] for n in range(var)]
    for j in range(0,var):
        gradii[j][j]=1
    def grad(ponto):
        grad=list()
        modular=list()
        tuple(ponto)
        for elemento in gradii:
            grad.append(diff(lambda x1,x2,x3,x4: g(x1,x2,x3,x4),ponto,tuple(elemento)))
        for element in grad:
            modular.append(element**2)
        modulo =  sum(modular)
        for i in range(0,var):
            grad[i]=grad[i]/sqrt( modulo)
        return grad
elif var == 5:
    gradii=[[0 for m in range(var)] for n in range(var)]
    for j in range(0,var):
        gradii[j][j]=1
    def grad(ponto):
        grad=list()
        modular=list()
        tuple(ponto)
        for elemento in gradii:
            grad.append(diff(lambda x1,x2,x3,x4,x5: g(x1,x2,x3,x4,x5),ponto,tuple(elemento)))
        for element in grad:
            modular.append(element**2)
        modulo =  sum(modular)
        for i in range(0,var):
            grad[i]=grad[i]/sqrt( modulo)
        return grad
## chute inicial: x
x=list()
x.append(raw_input('chute inicial: ' ) )
for line in x:
    line = line.split(",")
    x=map(float,line )
############################################
# definine o vetor step: (x - alpha.grad(x))
def vetor(alpha):
    zero=grad(x)
    mont=[0 for m in range(var)]
    for i in range(0,var):
        mont[i]=x[i] - alpha*zero[i]
    return mont
########### iteracao
k=0
steps=int(1000)
while k<= steps:
    alpha=[0,0,0,0]
#    print k
    if abs(max(grad(x),key=abs ))< epsilon:
        print "zero gradient, may have a minimum x:", x
        break
    else:
        a1=0
        a3=1
        g1=ng(x)
        g3=ng(vetor(a3))
        while ( g3 ) >= (g1):
            a3=a3*0.5
            g3= ng(vetor(a3) )
            if a3 < epsilon:
                print "nao tem mais o que fazer, talvez seja o minimo x:", x
                print "g(x):", g3
                break
        a2=a3*0.5
        g2= ng(vetor(a2))
# otmizacao do coeff alpha por polinomio de Newton
        h1=(g2-g1)/a2
        h2=(g3-g2)/(a3-a2)
        h3=(h2-h1)/a3
        a0=0.5*(a3 - (h1/h3))
        g0= ng(vetor(a0))
# encontrar o alpha para onde g(x) eh minimo
        gmin = [g0,g1,g2,g3]
        alpha = [a0,a1,a2,a3]
        ga = min(gmin)
        a = alpha[gmin.index(min(gmin))]
        if abs(ga - g1) < epsilon:
            mp.dps = 5; mp.pretty = True
            print "FEITO, x=", x
            print "steps", k
            print "g(x):", ga
            break
        x=vetor(a)
        print x
    k=k+1
    if k== steps:
        print "maximo de iteracoes atingido:", k
################################   FIM   #####################################
# example: (x1**3 +(x1**2)*x2 -x1*x3+6)**2 +(exp(x1)+exp(x2)-x3)**2 +(x2**2 -2*x1*x3 -4)**2
# -1.45536097438, -1.66625330383, 0.4204886158150
