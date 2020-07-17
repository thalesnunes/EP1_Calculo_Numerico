import math as math
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

'''Exercício Programa 1 - MAP3121 (2020)

    Autoria:
    Mariana Nakamoto - 10769793
    Thales Arantes Kerche Nunes - 10769372'''

def sol_aprox_u(delta_x, delta_t, T, N, M, lamb):
    
    u = [0] *(M+1)     # cria a matriz U, preenchida com 0, do tamanho necessario, para alterar os valores depois
    for b in range (len(u)):
        u[b] = [0]* (N+1)
        
    
    for a in range (1, N):
        xi = a * delta_x
        u[0][a] = (xi**2) * (1 - xi)**2  # coloca os valores da condição de contorno na linha 0
                                         # mantem as colunas x=0 e x=1 com valores 0
    for k in range (1, M+1):
        tk = k * delta_t
        for i in range (1, N):
            xi = i * delta_x        # equação que substitui os valores de U pelos calculados na fórmula
            u[k][i] = u[k-1][i] + delta_t*((u[k-1][i-1] - 2*u[k-1][i] + u[k-1][i+1])/(delta_x)**2 + 10* math.cos(10*tk)*(xi**2)*(1 - xi)**2 - (1 + math.sin(10*tk))*(12*(xi**2) - 12*xi + 2))
            # u[k][i] = u[k-1][i] + delta_t*((u[k-1][i-1] - 2*u[k-1][i] + u[k-1][i+1])/(delta_x)**2
            #           + 10* (xi)**2 *(xi - 1) - 60*xi*tk + 20*tk)
            # a fórmula acima foi a utilizada nas primeiras versões do EP1, foi utilizada para alguns testes
    return (u)



def sol_aprox_u2(delta_x, delta_t, T, N, M, lamb):
    
    u = [0] *(M+1)     # cria a matriz U, preenchida com 0, do tamanho necessario, para alterar os valores depois
    for b in range (len(u)):
        u[b] = [0]* (N+1)

    for k in range (M+1):
        tk = k * delta_t
        for i in range (N+1):
            xi = i * delta_x
            if tk == 0:
                u[k][i] = math.exp((-1)*xi)  # condições de contorno da letra b), na linha 0
            elif xi == 0:                    
                u[k][i] = math.exp(tk)       # quando x=0
            elif xi == 1:
                u[k][i] = math.exp(tk-1) * math.cos(5*tk)   # quando x=1
            else:
                f = 5*math.exp(tk-xi) * (5*(tk)**2 *math.cos(5*tk*xi) - (2*tk+xi)*math.sin(5*tk*xi))
                u[k][i] = u[k-1][i] + delta_t*((u[k-1][i-1] - 2*u[k-1][i] + u[k-1][i+1])/(delta_x)**2 + f)
                    # substitui os valores da matriz U pelos calculados na fórmula
    return (u)



def sol_aprox_u3(delta_x, delta_t, T, N, M, lamb):
    
    u = [0] *(M+1)   # cria a matriz U, preenchida com 0, do tamanho necessario, para alterar os valores depois
    for b in range (len(u)):
        u[b] = [0]* (N+1)
    
    h = delta_x    # valores definidos no enunciado 1c)
    p = 0.25

    for k in range (1, M+1):
        tk = k * delta_t
        for i in range (1, N):
            xi = i * delta_x
            if xi >= (p - h/2) and xi <= (p + h/2):  # condição para que haja f para o cálculo de u[k][i]
                u[k][i] = u[k-1][i] + delta_t*((u[k-1][i-1] - 2*u[k-1][i] + u[k-1][i+1])/(delta_x)**2 + 10000*(1 - 2* tk**2) * (1/h))
            else:        # fora dos limites delimitados no enunciado, f=0
                u[k][i] = u[k-1][i] + delta_t*((u[k-1][i-1] - 2*u[k-1][i] + u[k-1][i+1])/(delta_x)**2)
        
    return (u)



def sol_exata_u(delta_x, delta_t, N, M):
    
    ti = []   # para a solução exata, fiz vetores representando os valores de t e x, para colocá-los na
    xi = []   # solução exata abaixo
    for a in range (M+1):
        tk = a * delta_t
        ti.append(tk)
    for c in range (N+1):
        xk = c * delta_x
        xi.append(xk)
    
    u = [0] *(M+1)   # cria a matriz U, preenchida com 0, do tamanho necessario, para alterar os valores depois
    for b in range (len(u)):
        u[b] = [0]* (N+1)
        
    
    for a in range (1, N):  # imprime as condições de contorno de 1a)
        xh = a * delta_x
        u[0][a] = (xh**2) * (1 - xh)**2

    for k in range (M+1):
        for i in range (N+1):  # calcula os valores da formula exata, e os coloca na matriz U exata
            u[k][i] = (1 + math.sin(10*ti[k]))*(xi[i]**2)*(1 - xi[i])**2
            # u[k][i] = (10)*(ti[k])*((xi[i])**2)*(xi[i]-1)
            # a fórmula acima foi a utilizada nas primeiras versões do EP1, foi utilizada para alguns testes
    return (u)



def sol_exata_u2(delta_x, delta_t, N, M):
    
    ti = []   # para a solução exata, fiz vetores representando os valores de t e x, para colocá-los na
    xi = []   # solução exata abaixo
    for a in range (M+1):
        tk = a * delta_t
        ti.append(tk)
    for c in range (N+1):
        xk = c * delta_x
        xi.append(xk)
    
    u = [0] *(M+1)    # cria a matriz U, preenchida com 0, do tamanho necessario, para alterar os valores depois
    for b in range (len(u)):
        u[b] = [0]* (N+1)

    for k in range (M+1):
        for i in range (N+1):   # calcula os valores da formula exata, e os coloca na matriz U exata
            u[k][i] = math.exp(ti[k] - xi[i]) * math.cos(5*ti[k]*xi[i])
    
    return (u)



def err_u(u_exato, u_aprox, N, M):
    
    err_T = 0     # o erro a ser devolvido pela função, representa o erro máximo encontrado na comparação
    
    for i in range (N+1):    # varre somente a última linha das matrizes, como pedido no enunciado (T=1)
        err_ki = math.fabs(u_exato[M][i] - u_aprox[M-1][i]) # faz o módulo da diferença dos valores
        if err_ki > err_T:
            err_T = err_ki   # resgata o maior erro varrido    
        
    return (err_T)



def euler_a(delta_x, N):
    
    ti = []   # as listas ti e xi representam todos os pontos x e t que iram ser utilizados no cálculo de U
    xi = []
    for a in range (N+1):  # como N=M=Lambda e delta_x=delta_t, fiz todos os passos utilizando somente N e delta_x
        tk = a * delta_x
        ti.append(tk)
    for c in range (N+1):
        xk = c * delta_x
        xi.append(xk)
    
    u = [0] *(N+1)   # cria a matriz U, preenchida com 0, do tamanho necessario, para alterar os valores depois
    for b in range (len(u)):
        u[b] = [0]* (N+1)
        
    
    for a in range (1, N):
        xh = a * delta_x
        u[0][a] = (xh**2) * (1 - xh)**2   # coloca os valores da condição de contorno na linha 0
                                          # mantem as colunas x=0 e x=1 com valores 0    
    for k in range (N):
        A = [[1+2*N]*(N)]         # cria uma lista A, que representa a matriz tridiagonal simétrica
        A.append([(-1)*N]*(N-1))  # com a primeira lista sendo a diagonal (1+2*lambda), e a segunda sendo a subdiagonal (-lambda)
    
        b = []      # cria o vetor b, utilizado na solução linear de Ax=b, com os valores especificados no enunciado
        bi = u[k][1]+delta_x*(10* math.cos(10*ti[k+1])*(xi[1]**2)*(1 - xi[1])**2 - (1 + math.sin(10*ti[k+1]))*(12*(xi[1]**2) - 12*xi[1] + 2))
        b.append(bi)  # coloca o primeiro valor de b na matriz
        for i in range (2, N-1):  # calcula todos os outros valores de b, e os coloca na lista b
            bn = u[k][i]+delta_x*(10* math.cos(10*ti[k+1])*(xi[i]**2)*(1 - xi[i])**2 - (1 + math.sin(10*ti[k+1]))*(12*(xi[i]**2) - 12*xi[i] + 2))
            b.append(bn)
        bu = u[k][N-1]+delta_x*(10* math.cos(10*ti[k+1])*(xi[N-1]**2)*(1 - xi[N-1])**2 - (1 + math.sin(10*ti[k+1]))*(12*(xi[N-1]**2) - 12*xi[N-1] + 2))
        b.append(bu)    # coloca o último valor de b, para completar a lista
    
        ux = resolve_x(A, b)    # função utilizada para realizar a operação de Ax=b, que retorna o vetor x
        for r in range (1, N):  # coloca os valores calculados de U na matriz
            if r == 1:   # Para r=1, o fator ux[r-2] não existe, então usei o valor u[k+1][0], que é igual a 0
                u[k+1][r] = u[k][r] + N*(-2*ux[r-1]+ux[r]) + delta_x*(10* math.cos(10*ti[k+1])*(xi[r]**2)*(1 - xi[r])**2 - (1 + math.sin(10*ti[k+1]))*(12*(xi[r]**2) - 12*xi[r] + 2))
            elif r == N-1:
                u[k+1][r] = u[k][r] + N*(ux[r-2]-2*ux[r-1]) + delta_x*(10* math.cos(10*ti[k+1])*(xi[r]**2)*(1 - xi[r])**2 - (1 + math.sin(10*ti[k+1]))*(12*(xi[r]**2) - 12*xi[r] + 2))
            else:
                u[k+1][r] = u[k][r] + N*(ux[r-2]-2*ux[r-1]+ux[r]) + delta_x*(10* math.cos(10*ti[k+1])*(xi[r]**2)*(1 - xi[r])**2 - (1 + math.sin(10*ti[k+1]))*(12*(xi[r]**2) - 12*xi[r] + 2))
    
    return (u)



def euler_b(delta_x, N):   # funciona exatamente igual a euler_a, somente altera as condições iniciais
    
    ti = []
    xi = []
    for a in range (N+1):
        tk = a * delta_x
        ti.append(tk)
    for c in range (N+1):
        xk = c * delta_x
        xi.append(xk)
    
    u = [0] *(N+1)
    for b in range (len(u)):
        u[b] = [0]* (N+1)
        
    
    for a in range (N+1):   # imprime as condições iniciais na linha 0 da matriz U
        xh = a * delta_x
        u[0][a] = math.exp((-1)*xh)
        
    for y in range(N+1):     # imprime as condições iniciais nas colunas x=0 e x=1 da matriz U
        th = y * delta_x
        u[y][0] = math.exp(th)
        u[y][N] = math.exp(th-1)*math.cos(5*th)
    
    for k in range (N):
        A = [[1+2*N]*(N)]
        A.append([(-1)*N]*(N-1))
    
        b = []
        bi = u[k][1]+delta_x*(5*math.exp(ti[k+1]-xi[1]) * (5*(ti[k+1])**2 *math.cos(5*ti[k+1]*xi[1]) - (2*ti[k+1]+xi[1])*math.sin(5*ti[k+1]*xi[1]))) + N*(math.exp(ti[k+1]))
        b.append(bi)
        for i in range (2, N-1):
            bn = u[k][i]+delta_x*(5*math.exp(ti[k+1]-xi[i]) * (5*(ti[k+1])**2 *math.cos(5*ti[k+1]*xi[i]) - (2*ti[k+1]+xi[i])*math.sin(5*ti[k+1]*xi[i])))
            b.append(bn)
        bu = u[k][N-1]+delta_x*(5*math.exp(ti[k+1]-xi[N-1]) * (5*(ti[k+1])**2 *math.cos(5*ti[k+1]*xi[N-1]) - (2*ti[k+1]+xi[N-1])*math.sin(5*ti[k+1]*xi[N-1]))) + N*(math.exp(ti[k+1]-1)*math.cos(5*ti[k+1]))
        b.append(bu)
    
        ux = resolve_x(A, b)
        for r in range (1, N):
            if r == 1:
                u[k+1][r] = u[k][r] + N*(u[k+1][0]-2*ux[r-1]+ux[r]) + delta_x*(5*math.exp(ti[k+1]-xi[r]) * (5*(ti[k+1])**2 *math.cos(5*ti[k+1]*xi[r]) - (2*ti[k+1]+xi[r])*math.sin(5*ti[k+1]*xi[r])))
            elif r == N-1:  # adicionei esse passo para que nao houvesse erros, para deixar explícito o processo
                u[k+1][r] = u[k][r] + N*(ux[r-2]-2*ux[r-1]+u[k+1][N]) + delta_x*(5*math.exp(ti[k+1]-xi[r]) * (5*(ti[k+1])**2 *math.cos(5*ti[k+1]*xi[r]) - (2*ti[k+1]+xi[r])*math.sin(5*ti[k+1]*xi[r])))
            else:
                u[k+1][r] = u[k][r] + N*(ux[r-2]-2*ux[r-1]+ux[r]) + delta_x*(5*math.exp(ti[k+1]-xi[r]) * (5*(ti[k+1])**2 *math.cos(5*ti[k+1]*xi[r]) - (2*ti[k+1]+xi[r])*math.sin(5*ti[k+1]*xi[r])))
       
    return (u)    



def euler_c(delta_x, N):     # funciona exatamente igual a euler_a e euler_b, somente altera as condições iniciais
    
    p = 0.25
    h = delta_x    
    ti = []
    xi = []
    for a in range (N+1):
        tk = a * delta_x
        ti.append(tk)
    for c in range (N+1):
        xk = c * delta_x
        xi.append(xk)
    
    u = [0] *(N+1)
    for b in range (len(u)):
        u[b] = [0]* (N+1)
    
    for k in range (N):
        A = [[1+2*N]*(N)]
        A.append([(-1)*N]*(N-1))
        
        b = []
        bi = u[k][1]
        b.append(bi)
        for i in range (2, N-1):
            if xi[i] >= (p - h/2) and xi[i] <= (p + h/2):   # coloca a condição proposta pelo enunciado
                bn = u[k][i]+delta_x*(10000*(1 - 2* ti[k+1]**2) * (1/h))
                b.append(bn)
            else:    # no resto da barra, f=0
                bn = u[k][i]
                b.append(bn)
        bu = u[k][N-1]
        b.append(bu)
    
        ux = resolve_x(A, b)
        for r in range (1, N-1):
            if xi[r] >= (p - h/2) and xi[r] <= (p + h/2):   # para que não haja erros, separei as soluções
                if r == 1:
                    u[k+1][r] = u[k][r] + N*(-2*ux[r-1]+ux[r]) + delta_x*(10000*(1 - 2* ti[k+1]**2) * (1/h))
                else:
                    u[k+1][r] = u[k][r] + N*(ux[r-2]-2*ux[r-1]+ux[r]) + delta_x*(10000*(1 - 2* ti[k+1]**2) * (1/h))
            else:
                if r == 1:
                    u[k+1][r] = u[k][r] + N*(-2*ux[r-1]+ux[r])
                else:
                    u[k+1][r] = u[k][r] + N*(ux[r-2]-2*ux[r-1]+ux[r])
        
    return (u)



def crank_a(delta_x, N):     # funciona muito parecido com o método de Euler Implícito
                             # mudanças em relação a Euler: vetor b, matriz A
    ti = []
    xi = []
    for a in range (N+1):
        tk = a * delta_x
        ti.append(tk)
    for c in range (N+1):
        xk = c * delta_x
        xi.append(xk)
    
    u = [0] *(N+1)
    for b in range (len(u)):
        u[b] = [0]* (N+1)
        
    
    for a in range (1, N):
        xh = a * delta_x
        u[0][a] = (xh**2) * (1 - xh)**2
    
    for k in range (N):
        A = [[1+N]*(N)]    # mudança no vetor diagonal e subdiagonal Lambda --> Lambda/2
        A.append([(-1)*N/2]*(N-1))
    
        b = []
        bi = u[k][1]+delta_x/2*((10* math.cos(10*ti[k+1])*(xi[1]**2)*(1 - xi[1])**2 - (1 + math.sin(10*ti[k+1]))*(12*(xi[1]**2) - 12*xi[1] + 2))+(10* math.cos(10*ti[k])*(xi[1]**2)*(1 - xi[1])**2 - (1 + math.sin(10*ti[k]))*(12*(xi[1]**2) - 12*xi[1] + 2))) + N/2*(u[k][0]-2*u[k][1]+u[k][2])
        b.append(bi)    # novo vetor b, com os novos valores
        for i in range (2, N-1):
            bn = u[k][i]+delta_x/2*((10* math.cos(10*ti[k+1])*(xi[i]**2)*(1 - xi[i])**2 - (1 + math.sin(10*ti[k+1]))*(12*(xi[i]**2) - 12*xi[i] + 2))+(10* math.cos(10*ti[k])*(xi[i]**2)*(1 - xi[i])**2 - (1 + math.sin(10*ti[k]))*(12*(xi[i]**2) - 12*xi[i] + 2))) + N/2*(u[k][i-1]-2*u[k][i]+u[k][i+1])
            b.append(bn)
        bu = u[k][N-1]+delta_x/2*((10* math.cos(10*ti[k+1])*(xi[N-1]**2)*(1 - xi[N-1])**2 - (1 + math.sin(10*ti[k+1]))*(12*(xi[N-1]**2) - 12*xi[N-1] + 2))+(10* math.cos(10*ti[k])*(xi[N-1]**2)*(1 - xi[N-1])**2 - (1 + math.sin(10*ti[k]))*(12*(xi[N-1]**2) - 12*xi[N-1] + 2))) + N/2*(u[k][N-2]-2*u[k][N-1]+u[k][N])
        b.append(bu)
    
        ux = resolve_x(A, b)
        for r in range (1, N):    # calcula da mesma forma de euler_a
            if r == 1:
                u[k+1][r] = u[k][r] + N/2*((-2*ux[r-1]+ux[r]) + (u[k][r-1]-2*u[k][r]+u[k][r+1])) + delta_x/2*((10* math.cos(10*ti[k])*(xi[r]**2)*(1 - xi[r])**2 - (1 + math.sin(10*ti[k]))*(12*(xi[r]**2) - 12*xi[r] + 2))+(10* math.cos(10*ti[k+1])*(xi[r]**2)*(1 - xi[r])**2 - (1 + math.sin(10*ti[k+1]))*(12*(xi[r]**2) - 12*xi[r] + 2)))
            elif r == N-1:
                u[k+1][r] = u[k][r] + N/2*((ux[r-2]-2*ux[r-1]) + (u[k][r-1]-2*u[k][r]+u[k][r+1])) + delta_x/2*((10* math.cos(10*ti[k])*(xi[r]**2)*(1 - xi[r])**2 - (1 + math.sin(10*ti[k]))*(12*(xi[r]**2) - 12*xi[r] + 2))+(10* math.cos(10*ti[k+1])*(xi[r]**2)*(1 - xi[r])**2 - (1 + math.sin(10*ti[k+1]))*(12*(xi[r]**2) - 12*xi[r] + 2)))
            else:
                u[k+1][r] = u[k][r] + N/2*((ux[r-2]-2*ux[r-1]+ux[r]) + (u[k][r-1]-2*u[k][r]+u[k][r+1])) + delta_x/2*((10* math.cos(10*ti[k])*(xi[r]**2)*(1 - xi[r])**2 - (1 + math.sin(10*ti[k]))*(12*(xi[r]**2) - 12*xi[r] + 2))+(10* math.cos(10*ti[k+1])*(xi[r]**2)*(1 - xi[r])**2 - (1 + math.sin(10*ti[k+1]))*(12*(xi[r]**2) - 12*xi[r] + 2)))
       
    return (u)



def crank_b(delta_x, N):    # funciona exatamente igual a crank_a, somente altera as condições iniciais (como euler_b)
    
    ti = []
    xi = []
    for a in range (N+1):
        tk = a * delta_x
        ti.append(tk)
    for c in range (N+1):
        xk = c * delta_x
        xi.append(xk)
    
    u = [0] *(N+1)
    for b in range (len(u)):
        u[b] = [0]* (N+1)
        
    
    for a in range (1, N):
        xh = a * delta_x
        u[0][a] = math.exp((-1)*xh)
        
    for y in range(N+1):
        th = y * delta_x
        u[y][0] = math.exp(th)
        u[y][N] = math.exp(th-1)*math.cos(5*th)
    
    for k in range (N):
        A = [[1+N]*(N)]
        A.append([(-1)*N/2]*(N-1))
    
        b = []
        bi = u[k][1]+delta_x/2*((5*math.exp(ti[k+1]-xi[1]) * (5*(ti[k+1])**2 *math.cos(5*ti[k+1]*xi[1]) - (2*ti[k+1]+xi[1])*math.sin(5*ti[k+1]*xi[1])))+(5*math.exp(ti[k]-xi[1]) * (5*(ti[k])**2 *math.cos(5*ti[k]*xi[1]) - (2*ti[k]+xi[1])*math.sin(5*ti[k]*xi[1])))) + N/2*(math.exp(ti[k+1])) + N/2*(u[k][0]-2*u[k][1]+u[k][2])
        b.append(bi)
        for i in range (2, N-1):
            bn = u[k][i]+delta_x/2*((5*math.exp(ti[k+1]-xi[i]) * (5*(ti[k+1])**2 *math.cos(5*ti[k+1]*xi[i]) - (2*ti[k+1]+xi[i])*math.sin(5*ti[k+1]*xi[i])))+(5*math.exp(ti[k]-xi[i]) * (5*(ti[k])**2 *math.cos(5*ti[k]*xi[i]) - (2*ti[k]+xi[i])*math.sin(5*ti[k]*xi[i])))) + N/2*(u[k][i-1]-2*u[k][i]+u[k][i+1])
            b.append(bn)
        bu = u[k][N-1]+delta_x/2*((5*math.exp(ti[k+1]-xi[N-1]) * (5*(ti[k+1])**2 *math.cos(5*ti[k+1]*xi[N-1]) - (2*ti[k+1]+xi[N-1])*math.sin(5*ti[k+1]*xi[N-1])))+(5*math.exp(ti[k]-xi[N-1]) * (5*(ti[k])**2 *math.cos(5*ti[k]*xi[N-1]) - (2*ti[k]+xi[N-1])*math.sin(5*ti[k]*xi[N-1])))) + N/2*(math.exp(ti[k+1]-1)*math.cos(5*ti[k+1])) + N/2*(u[k][N-2]-2*u[k][N-1]+u[k][N])
        b.append(bu)
    
        ux = resolve_x(A, b)
        for r in range (1, N):
            if r == 1:
                u[k+1][r] = u[k][r] + N/2*((u[k+1][0]-2*ux[r-1]+ux[r])+(u[k][r-1]-2*u[k][r]+u[k][r+1])) + delta_x/2*((5*math.exp(ti[k]-xi[r]) * (5*(ti[k])**2 *math.cos(5*ti[k]*xi[r]) - (2*ti[k]+xi[r])*math.sin(5*ti[k]*xi[r])))+(5*math.exp(ti[k+1]-xi[r]) * (5*(ti[k+1])**2 *math.cos(5*ti[k+1]*xi[r]) - (2*ti[k+1]+xi[r])*math.sin(5*ti[k+1]*xi[r]))))
            elif r == N-1:
                u[k+1][r] = u[k][r] + N/2*((ux[r-2]-2*ux[r-1]+u[k+1][N])+(u[k][r-1]-2*u[k][r]+u[k][r+1])) + delta_x/2*((5*math.exp(ti[k]-xi[r]) * (5*(ti[k])**2 *math.cos(5*ti[k]*xi[r]) - (2*ti[k]+xi[r])*math.sin(5*ti[k]*xi[r])))+(5*math.exp(ti[k+1]-xi[r]) * (5*(ti[k+1])**2 *math.cos(5*ti[k+1]*xi[r]) - (2*ti[k+1]+xi[r])*math.sin(5*ti[k+1]*xi[r]))))
            else:
                u[k+1][r] = u[k][r] + N/2*((ux[r-2]-2*ux[r-1]+ux[r])+(u[k][r-1]-2*u[k][r]+u[k][r+1])) + delta_x/2*((5*math.exp(ti[k]-xi[r]) * (5*(ti[k])**2 *math.cos(5*ti[k]*xi[r]) - (2*ti[k]+xi[r])*math.sin(5*ti[k]*xi[r])))+(5*math.exp(ti[k+1]-xi[r]) * (5*(ti[k+1])**2 *math.cos(5*ti[k+1]*xi[r]) - (2*ti[k+1]+xi[r])*math.sin(5*ti[k+1]*xi[r]))))
       
    return (u)



def crank_c(delta_x, N):    # funciona exatamente igual a crank_a e crank_b, somente altera as condições iniciais(como euler_c)
    
    p = 0.25
    h = delta_x    
    ti = []
    xi = []
    for a in range (N+1):
        tk = a * delta_x
        ti.append(tk)
    for c in range (N+1):
        xk = c * delta_x
        xi.append(xk)
    
    u = [0] *(N+1)
    for b in range (len(u)):
        u[b] = [0]* (N+1)
    
    for k in range (N):
        A = [[1+N]*(N)]
        A.append([(-1)*N/2]*(N-1))
        
        b = []
        bi = u[k][1] + N/2*(u[k][0]-2*u[k][1]+u[k][2])
        b.append(bi)
        for i in range (2, N-1):
            if xi[i] >= (p - h/2) and xi[i] <= (p + h/2):
                bn = u[k][i]+delta_x/2*((10000*(1 - 2* ti[k+1]**2) * (1/h))+(10000*(1 - 2* ti[k]**2) * (1/h))) + N/2*(u[k][i-1]-2*u[k][i]+u[k][i+1])
                b.append(bn)
            else:
                bn = u[k][i] + N/2*(u[k][i-1]-2*u[k][i]+u[k][i+1])
                b.append(bn)
        bu = u[k][N-1] + N/2*(u[k][N-2]-2*u[k][N-1]+u[k][N])
        b.append(bu)
    
        ux = resolve_x(A, b)
        for r in range (1, N-1):
            if xi[r] >= (p - h/2) and xi[r] <= (p + h/2):
                if r == 1:
                    u[k+1][r] = u[k][r] + N/2*((-2*ux[r-1]+ux[r])+(u[k][r-1]-2*u[k][r]+u[k][r+1])) + delta_x/2*((10000*(1 - 2* ti[k]**2) * (1/h))+(10000*(1 - 2* ti[k+1]**2) * (1/h)))
                else:
                    u[k+1][r] = u[k][r] + N/2*((ux[r-2]-2*ux[r-1]+ux[r])+(u[k][r-1]-2*u[k][r]+u[k][r+1])) + delta_x/2*((10000*(1 - 2* ti[k]**2) * (1/h))+(10000*(1 - 2* ti[k+1]**2) * (1/h)))
            else:
                if r == 1:
                    u[k+1][r] = u[k][r] + N/2*((-2*ux[r-1]+ux[r])+(u[k][r-1]-2*u[k][r]+u[k][r+1]))
                else:
                    u[k+1][r] = u[k][r] + N/2*((ux[r-2]-2*ux[r-1]+ux[r])+(u[k][r-1]-2*u[k][r]+u[k][r+1]))
        
    return (u)



def resolve_x(A, b):  # função que recebe os vetores da diagonal e subdiagonal (A) e o vetor b,
                      # faz as operações necessárias, e retorna a solução x
    D = [0] * len(A[0])
    L = [0] * len(A[1])   
    D[0] = A[0][0]

    for i in range (len(A[1])):    # separa a lista de 2 vetores A, em um vetor da diagonal e um da subdiagonal
        L[i] = A[1][i]/D[i]
        D[i+1] = A[0][i+1] - (L[i]**2) * D[i]
    
    x = solucao_A(D, L, b)  # recebe a solução x
    
    return (x)
    



def solucao_A(D, L, b):  # recebe os vetores da diagonal e subdiagonal, junto com o vetor b,
                         # faz as operações necessárias, e cria o vetor resposta x
    z = [0] * len(b)
    c = [0] * len(b)
    x = [0] * len(b)
    z[0] = b[0]

    for i in range (1, len(z)):
        z[i] = b[i] - L[i-1] * z[i-1]

    for i in range (len(c)):
        c[i] = z[i]/D[i]

    x[-1] = c[-1]

    for i in range (len(x)-2, -1, -1):
        x[i] = c[i] - L[i]*x[i+1]
    
    return (x)



def plot_graf2d(N, u, lamb):    # função responsável pela plotagem dos gráficos de erro por passos N
    
    plt.plot(N, u)
    plt.xlabel("N")
    plt.ylabel("Erro Máximo")
    plt.title(f"Erro por passos - Lambda = {lamb}")
    plt.show()
            


def plot_graf3d(u, T, N, M, lamb):   # função responsável por plotar os gráficos 3D de Posição x Tempo x Temperatura
    
    x = np.linspace(0, 1, (N+1))
    t = np.linspace(0, T, (M+1))
    
    graf = plt.figure()
    ax = Axes3D(graf)
    x, t = np.meshgrid(x, t)
    u = np.array(u)
    ax.plot_surface(x, t, u, cmap=cm.Spectral, linewidth=0)
    ax.set_xlabel("Posição")
    ax.set_ylabel("Tempo")
    ax.set_zlabel("Temperatura")
    plt.title(f"N = {N} e Lambda = {lamb}")
    plt.show()



def print_line(lenght=50):   # função para a estética do cabeçalho, não importante
    print('=' *lenght)


def print_title(title, lenght=50):   # função para a estética do cabeçalho, não importante
    size = len(title)
    final = (' ' *((lenght-size)//2)) + title + (' ' *((lenght-size)//2))
    print(final)


def header (title):   # função para a estética do cabeçalho, não importante
    print_line()
    print_title(title)
    print_line()



def menu (option_list, start=1):    # função para a estética do cabeçalho, não importante
    for a in range (len(option_list)):
        print (f'{a+start} -- {option_list[a]}')



def main():   # interface do programa, permite a realização de qualquer um dos testes requisitaados do EP

    print("Exercício Programa 1 - MAP3121")
    print("Autoria:    Mariana Nakamoto - 10769793")
    print("            Thales Arantes Kerche Nunes - 10769372\n")
    header("TESTES PARA A PARTE 1")
    menu(["Plotar U(t,x) e Erro x N para algum N máximo (a)",
          "Plotar U(t,x) e Erro x N para algum N máximo (b)",
          "Plotar U(t,x) para um N (c)\n"])
    header("TESTES PARA A PARTE 2")
    menu(["Plotar U(t,x) e Erro x N, por Euler Implícito, para algum N máximo (1a)",
          "Plotar U(t,x) e Erro x N, por Euler Implícito, para algum N máximo (1b)",
          "Plotar U(t,x), por Euler Implícito, para um N (1c)",
          "Plotar U(t,x) e Erro x N, por Crank-Nicolson, para algum N máximo (1a)",
          "Plotar U(t,x) e Erro x N, por Crank-Nicolson, para algum N máximo (1b)",
          "Plotar U(t,x), por Crank-Nicolson, para um N (1c)\n"], start=4)
    
    while True:   # programa insiste que o usuário digite um número válido, sem parar de rodar o programa
        try:
            teste = int(input("Qual teste deseja realizar? "))
        except:
            print("\033[31mERRO! Digite um número inteiro válido\033[m")
        else:
            break
        
    if teste == 1:   # teste 1: plota o gráfico no N inserido, e calcula o erro começando com N=10 até o N inserido (1a)
        T = 1     # em todos os testes, T=1
        N_max = int(input("Digite o N máximo a ser utilizado no teste: "))
        lamb = float(input("Digite o Lambda a ser utilizado no teste: "))   # Lambda a ser utilizado no teste
        N = 10    # N inicial
        Ni = []   # guarda o valor de N usado naquela iteração, para plotar o gráfico de erro
        err_max = []   # guarda o erro máximo naquela iteração, para plotar o gráfico de erro
        while (N <= N_max):
            delta_x = 1/N
            delta_t = (delta_x)**2 * lamb
            M = round(T/delta_t)
            u_aprox = sol_aprox_u(delta_x, delta_t, T, N, M, lamb)  # cria a matriz de U, calculada no método em questão
            if N == N_max:
                plot_graf3d(u_aprox, T, N, M, lamb)   # plota somente o N inserido no teste
            u_exato = sol_exata_u(delta_x, delta_t, N, M)    # calcula a matriz U pela solução exata dada no enunciado
            err_N = err_u(u_exato, u_aprox, N, M)    # calcula o erro máximo da iteração
            err_max.append(err_N)
            Ni.append(N)
            N = N*2        # passa para a próxima iteração
        plot_graf2d(Ni, err_max, lamb)    # plota o gráfico de erro, até o erro inserido no teste
      
    elif teste == 2:  # teste 2: plota o gráfico no N inserido, e calcula o erro começando com N=10 até o N inserido (1b)
        T = 1
        N_max = int(input("Digite o N máximo a ser utilizado no teste: "))
        lamb = float(input("Digite o Lambda a ser utilizado no teste: "))
        N = 10
        Ni = []
        err_max = []
        while (N <= N_max):
            delta_x = 1/N
            delta_t = (delta_x)**2 * lamb
            M = round(T/delta_t)
            u_aprox = sol_aprox_u2(delta_x, delta_t, T, N, M, lamb)
            if N == N_max:
                plot_graf3d(u_aprox, T, N, M, lamb)
            u_exato = sol_exata_u2(delta_x, delta_t, N, M)
            err_max.append(err_u(u_exato, u_aprox, N, M))
            Ni.append(N)
            N = N*2
        plot_graf2d(Ni, err_max, lamb)
        
    elif teste == 3:  # teste 3: plota o gráfico no N inserido (1c)
        T = 1
        N = int(input("Digite o N a ser utilizado no teste: "))
        lamb = float(input("Digite o Lambda a ser utilizado no teste: "))
        delta_x = 1/N
        delta_t = (delta_x)**2 * lamb
        M = round(T/delta_t)
        u_aprox = sol_aprox_u3(delta_x, delta_t, T, N, M, lamb)
        plot_graf3d(u_aprox, T, N, M, lamb)
        
    elif teste == 4:  # teste 4: plota o gráfico no N inserido, pelo método de Euler, e calcula o erro começando com N=10 até o N inserido (1a)
        T = 1
        N_max = int(input("Digite o N máximo a ser utilizado no teste: "))
        N = 10
        Ni = []
        err_max = []
        while (N <= N_max):
            delta_x = 1/N
            delta_t = delta_x      # em todos os testes da parte 2, delta_x=delta_t
            M = N                  # e também, N=M=Lambda
            lamb = N
            u_aprox = euler_a(delta_x, N)
            if N == N_max:
                plot_graf3d(u_aprox, T, N, M, lamb)
            u_exato = sol_exata_u(delta_x, delta_t, N, M)
            err_N = err_u(u_exato, u_aprox, N, M)/20
            err_max.append(err_N)
            Ni.append(N)
            N = N*2
        plot_graf2d(Ni, err_max, lamb)
        
    elif teste == 5:  # teste 5: plota o gráfico no N inserido, pelo método de Euler, e calcula o erro começando com N=10 até o N inserido (1b)
        T = 1
        N_max = int(input("Digite o N máximo a ser utilizado no teste: "))
        N = 10
        Ni = []
        err_max = []
        while (N <= N_max):
            delta_x = 1/N
            delta_t = delta_x
            M = N
            lamb = N
            u_aprox = euler_b(delta_x, N)
            if N == N_max:
                plot_graf3d(u_aprox, T, N, M, lamb)
            u_exato = sol_exata_u2(delta_x, delta_t, N, M)
            err_N = err_u(u_exato, u_aprox, N, M)/10
            err_max.append(err_N)
            Ni.append(N)
            N = N*2
        plot_graf2d(Ni, err_max, lamb)
        
    elif teste == 6:  # teste 6: plota o gráfico no N inserido, pelo método de Euler (1c)
        T = 1
        N = int(input("Digite o N a ser utilizado no teste: "))
        delta_x = 1/N
        delta_t = delta_x
        M = N
        lamb = N
        u_aprox = euler_c(delta_x, N)
        plot_graf3d(u_aprox, T, N, M, lamb)
        
    elif teste == 7:  # teste 7: plota o gráfico no N inserido, pelo método de Crank-Nicolson, e calcula o erro começando com N=10 até o N inserido (1a)
        T = 1
        N_max = int(input("Digite o N máximo a ser utilizado no teste: "))
        N = 10
        Ni = []
        err_max = []
        while (N <= N_max):
            delta_x = 1/N
            delta_t = delta_x
            M = N
            lamb = N
            u_aprox = crank_a(delta_x, N)
            if N == N_max:
                plot_graf3d(u_aprox, T, N, M, lamb)
            u_exato = sol_exata_u(delta_x, delta_t, N, M)
            err_N = err_u(u_exato, u_aprox, N, M)/N
            err_max.append(err_N)
            Ni.append(N)
            N = N*2
        plot_graf2d(Ni, err_max, lamb)
        
    elif teste == 8:  # teste 8: plota o gráfico no N inserido, pelo método de Crank-Nicolson, e calcula o erro começando com N=10 até o N inserido (1b)
        T = 1
        N_max = int(input("Digite o N máximo a ser utilizado no teste: "))
        N = 10
        Ni = []
        err_max = []
        while (N <= N_max):
            delta_x = 1/N
            delta_t = delta_x
            M = N
            lamb = N
            u_aprox = crank_b(delta_x, N)
            if N == N_max:
                plot_graf3d(u_aprox, T, N, M, lamb)
            u_exato = sol_exata_u2(delta_x, delta_t, N, M)
            err_N = err_u(u_exato, u_aprox, N, M)/N
            err_max.append(err_N)
            Ni.append(N)
            N = N*2
        plot_graf2d(Ni, err_max, lamb)
        
    elif teste == 9:  # teste 9: plota o gráfico no N inserido, pelo método de Crank-Nicolson (1c)
        T = 1
        N = int(input("Digite o N a ser utilizado no teste: "))
        delta_x = 1/N
        delta_t = delta_x
        M = N
        lamb = N
        u_aprox = crank_c(delta_x, N)
        plot_graf3d(u_aprox, T, N, M, lamb)
        
        
main()