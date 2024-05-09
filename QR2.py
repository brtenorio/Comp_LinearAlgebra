#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import sys

######################################################################
# Metodo de Householder com QR,
# para calculo de autovalores de matrizes simetricas e nao simetricas
# Bruno Nunes Cabral Tenorio
# Maio de 2016
######################################################################
#

#
def input_matrix():
    '''
    Funcao usada para entrada da matrix pelo usuario
    output:     Matriz
    '''
        
    din= int(raw_input('dimensao da matriz: ' ))
    matrix = list()
    
    # coleta dos dados do sistema
    linhas = list()
    for n in range(0,din):
        linhas.append(raw_input('linha %s da matriz: ' %(n +1) ))
    for line in linhas:
        line = line.split(",")
        matrix.append(map(float,line ))
    
    return np.matrix(matrix)

def modulo(vetor):
    '''
    Funcao que calcula o modulo
    input:      Vetor
    output:     Valor do modulo em float
    '''
    
    return np.sqrt(np.sum(np.multiply(vetor, vetor)))

def Householder(matrix, hholder):
    '''
    Funcao que realiza o processo de iteracao de Householder
    input:      matrix -> Matriz do problema      
                hholder -> Matriz upper-Hessemberg ou tridiagonal
    output:     matrix -> Matriz do problema transformada
    '''
       
    din = len(matrix)
    ################################### proc. de iteracao Householder
    for k in range(0, din-2):
        if hholder == 'sim':
            pass
        elif hholder == 'nao':
    # matriz householder=Iholder
            Iholder = np.zeros((din,din))
            for i in range(0,din):
                Iholder[i,i]=1.
    # matriz reduzida A(k):
            ak = np.zeros((din-1-k, din-1-k))
            for i in range(0,din-1-k):
                for j in range(0,din-1-k):
                    ak[i,j]=matrix[i+1+k,j+1+k]
    # definicoes dos vetores e matrizes
            u = np.zeros((din-1-k))
            v = np.zeros((din-1-k))
            w = np.zeros((din-1-k))
            ident = np.identity(din-1-k)
            wwt = np.zeros((din-1-k,din-1-k))
    # vetor u e v
            for i in xrange(0,din-1-k):
                u[i]=matrix[i+k+1,k]
            v[0] = modulo(u)
    # vetor w
            for i in xrange(0,din-1-k):
                w[i]=u[i]-v[i]
            modW=modulo(w)
            for i in xrange(0,din-1-k):
                w[i]=w[i]/modW
    # matriz wwt
            for i in xrange(0,din-1-k):
                for j in xrange(0,din-1-k):
                    wwt[i,j]=w[i]*w[j]
    # matrix de reflexao de holder
            holder = np.add(ident, -2*wwt)
    # Iholder -> matriz householder de dimencao completa
    # lida tambem com matrizes nao simetricas
            for m in range(0,din-1-k):
                for n in range(0,din-1-k):
                    Iholder[m+1+k][n+1+k]=holder[m][n]
    # produto de matrizes; Iholder.matrx.Iholder
            AP = np.dot(matrix,Iholder)
            matrix = np.dot(Iholder,AP)
    ############ END

    # elimina valores muito pequenos de matrix
    for i in xrange(0,din):
        for j in xrange(0,din):
            matrix[i,j]=round(matrix[i,j],8)
            if np.abs(matrix[i,j]) < 10**-8:
                matrix[i,j] = 0.
#    for rq in matrix:
#        print("matriz upper-Hessemberg Linha %s :"  %(matrix.index(rq)+1),  rq)
#    print("fim do metodo de Householder")
    return matrix
#
def QR(matrix):
    '''
    Funcao que realiza o processo QR na matrix dada
    input:  Matrix -> Matriz que sofrera o processo
    output: 
    '''
    hholder= raw_input('a matriz eh upper-Hessemberg (ou tridiagonal)?: sim ou nao ' )
    din = len(matrix)
    matrix = Householder(matrix, hholder)
    print("matriz upper-Hessemberg (ou tridiagonal)")
    print(matrix)
    
    ############## PROCESSO QR ##########################################
    ################################### proc. de iteracao
    # limites do proc. : etpas= 200, de iteracao=10**-4
    k=0
    shift =0
    while k< 200:
    # Qn vai acumular as matrizes de rotacao Pn
        Qn=list()
    # teste de parada do loop -> limite = 10**-4
        lista_parada=list()
        for i in range(0,din-1):
            lista_parada.append(matrix[i+1,i])
        if abs(max(lista_parada,key=abs)) < 10**-4:
            break
    ## calculo do shift
        for a in xrange(0,din-1):
            if matrix[a+1,a] > 0.01:
                bs= (matrix[a+1,a+1] + matrix[a,a])
                cs=matrix[a+1,a+1]*matrix[a,a]-matrix[a+1,a]*matrix[a,a+1]
                ds= bs**2 - 4*cs
                if ds >= 0:
                    mi1=(bs+ (ds)**0.5)/2
                    mi2=(bs- (ds)**0.5)/2
                else:
                   continue
                lambdalista=np.array((mi1,mi2))
                lambdamin=np.array((mi1 - np.abs(matrix[a,a]), mi2-np.abs(matrix[a,a])))
                theta=lambdalista[np.argmin(lambdamin)]
                shift = shift + theta
                for m in range(0,din):
                    matrix[m,m] = matrix[m,m]- theta
                break 
    ## fim do shift
    ### pegando os valores a,b e construindo as matrizes P e R
        for j in xrange(0, din-1):
    # matriz P
            Pn=np.identity(din)
    #
            if np.abs(matrix[j+1,j]) < 10**-7:
                continue
    # verifica se b(k+1) ja eh pequeno e constroi P 
            elif np.abs(matrix[j+1,j]) > 10**-7:
                s= matrix[j+1,j]/(((matrix[j+1,j])**2+(matrix[j,j])**2)**0.5)
                c= matrix[j,j]/(((matrix[j+1,j])**2+(matrix[j,j])**2)**0.5)
                Pn[j,j]=Pn[j+1,j+1]=c
                Pn[j+1,j]=-s
                Pn[j,j+1]=s
    # R(n) = P(n).A(n)
                R=np.dot(Pn,matrix)
    # eliminando de R valores muito pequenos
                for x in xrange(0,din):
                    for y in xrange(0,din):
                        if np.abs(R[x,y]) <= 10**-9:
                            R[x,y] = 0.
                matrix = R
    # ja estou armazenando em Qn os Pn(transposto)
                Qn.append(zip(*Pn))
    # Qt=Qn[0].Qn[1]...
        Qt=Qn[0]
        for v in xrange(1,len(Qn)):
            Qt=np.dot(Qt,Qn[v])            
    # A(n+1)=R(n).Q(n)
        matrix=np.dot(R,Qt)
        k=k+1
    #    print k
    #
    # eliminando valores muito pequenos de A
    for x in xrange(0,din):
        for y in xrange(0,din):
            if np.abs(matrix[x,y]) <= 10**-6:
                matrix[x,y] = 0.
    
    ####### imprime o resultado
    print("Fim do metodo QR")
    print( matrix)
#    for rq in matrix:
#        print("matriz quasi-Diagonal Linha %s :"  %(matrix.index(rq)+1),  rq)
#    #
    for i in range(0,din):
        print("autovalor", matrix[i,i] + shift)
    print("QR steps:", k)
    #
    ################################# END
    
#### ESSA ALTERNATIVA IMPRIME EM UM ARQUIVO DE TEXTO "output.txt"
##with open('output.txt', 'w') as f:
##    sys.stdout = f
##    print("Fim do metodo QR")
##    for rq in matrix:
##        print("matriz quasi-Diagonal Linha %s :"  %(matrix.index(rq)+1),  rq)
####
##    for i in range(0,din):
##        print("autovalor", matrix[i][i])
##    print("QR steps:", k)

matrix = input_matrix()
QR(matrix)






