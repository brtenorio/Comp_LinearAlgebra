######################################################################
# Metodo de Gauss-Jacobi
# Bruno Nunes Cabral Tenorio
# Abril de 2016
######################################################################

din= raw_input('dimensao do sistema linear: ' )
din = int(din)

epson=raw_input('tolerancia erro: ' )
epson=float(epson)

Nint=raw_input('N maximo de iteracoes: ' )
Nint=int(Nint)

# chute inicial
X=list()
X.append(raw_input('chute inicial: ' ) )
for line in X:
    line = line.split(",")
    X=map(float,line )
#

matrix=list()
b=list()
# coleta dos dados do sistema
# matrix
linhas=list()
for n in range(0,din):
    linhas.append(raw_input('linha %s do sistema linear: ' %(n +1) ))
for line in linhas:
    line = line.split(",")
    matrix.append(map(float,line ))
# vetor b
solucao=list()
solucao.append(raw_input('vetor solucao escrito como linha: ' ) )
for line in solucao:
    line = line.split(",")
    b=map(float,line )
# fim da coleta de dados

print 'matriz normal: '
print matrix
print 'vetor solucao (normal) transposto: '
print b

###pivoteamento parcial
### faz o pivoteamento na coluna v
##def pivo(v):
##    vetor=list()
##    for i in range(v,len(matrix)):
##        vetor.append(matrix[i][v])
##    ind=vetor.index(max(vetor, key=abs))
##    matrix[v],matrix[ind+v] = matrix[ind+v],matrix[v]
##    b[v],b[ind+v]=b[ind+v],b[v]
##    return matrix,b
#####
##
### metodo de eliminacao de Gauss 
##m = [[0 for x in range(din)] for y in range(din)]
##
##for j in range(0, din-1):
###    print 'teste', j
##    pivo(j)
###    print matrix, b
##    for k in range(int(1+j), din ):
###        print k
##        m[j][k-j-1] = matrix[k][j]/matrix[j][j]
##        b[k]=b[k]-m[j][k-j-1]*b[j]
##        for i in  range(j, din):
##            matrix[k][i] = matrix[k][i] - m[j][k-j-1]*matrix[j][i]
##    for r in range(0,len(matrix)):
##        for t in range(0,din ):
##            if abs( matrix[r][t] ) < 10**-8:
##                matrix[r][t]= 0.0
##        if abs( b[r] ) < 10**-8:
##            b[r] = 0.0
##
###            a.append( float(matrix[k][i]) - m[k-1]*float(matrix[j][i]) )
###    print m[j]
##
##print 'matriz triangular superior: '
##print matrix
##print 'vetor solucao transposto: '
##print b
##print 'multiplicadores: '
##print m
##
### Solucoes x_i do sistema
##x=[0 for x in range(din)]
##x[din-1]=b[din-1]/matrix[din-1][din-1]
##
##def aux(p):
##    ax=list()
##    for q in range(din-1,p-1,-1):
##        ax.append(matrix[p][q]*x[q])
##    return ax
##for h in range(din-2,-1,-1):
##    x[h]=(b[h]-sum(aux(h) ) )/matrix[h][h]
##print 'solucoes do sistema linear, vetor X'
##for xi in x:
##    print "x _ %s :" %(x.index(xi)+1), round(xi,3)
#########################################################################

##  Metodo iterativo de Gauss Jacobi

# Separar a matrix em D-L-U

D = [[0 for x in range(din)] for y in range(din)]

L = [[0 for x in range(din)] for y in range(din)]

U = [[0 for x in range(din)] for y in range(din)]

C = list()

#
for i in range( 0, din ):
    D[i][i]=matrix[i][i]
print "D =", D
#
for i in range(0, din ):
    for j in range( 0, din ):
        if i>j :
            L[i][j]=(-1)*matrix[i][j]

print "L =", L
#
for i in range(0, din ):
    for j in range( 0, din ):
        if i<j :
            U[i][j]=(-1)*matrix[i][j]

print "U =", U
#
T = [[0 for x in range(din)] for y in range(din)]
for k in range(0, din):
    for l in range(0, din):
        T[k][l]= (1/D[k][k])*(L[k][l]+U[k][l])
print "T = ", T
#
for i in range(0, din):
    C.append((1/D[i][i])*b[i])
print "C =", C

## definicao de uma funcao (dot) que faz produto matriz vetor mais vetor
## => dot(A,B,G) = A.B + G
## A eh matrix, B e G vetores (G pode ser 0)
def dot(A,B,G):
    vec=list()
    aux=[[0 for x in range(din)] for y in range(din)]
    for i in range(0, din):
        for j in range(0, din):
            aux[i][j] = A[i][j]*B[j]
    for k in range( 0, din ):
        if G==0:
            vec.append( sum(aux[k]) )
        else:
            vec.append( sum(aux[k]) + G[k] )
    return vec

### iteracao do metodo
n=0
while n < Nint:
    norm=[0 for x in range(din)]
    ak=X
    X=dot(T,X,C)
    for i in range(0,din):
        norm[i] = ak[i]-X[i]
    if abs(max(norm, key=abs )) < epson:
        break
    else:
        n=n+1

print "solucao X =", X
print "erro, Nint :", max(norm, key=abs), n
