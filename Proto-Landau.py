######################################################################################################

# "Achador" Numérico de Autoenergias das Ondas de Spin

######################################################################################################

# Importação de Módulos

import numpy as np
import matplotlib.pyplot as plt

######################################################################################################

# Desenvolvimento

pi = np.pi

a0 = float(input('a0 = '))   # Escolhendo um valor para o parâmetro de rede;
J = float(input('J = '))   # Escolhendo a constante de exchange;
S = float(input('S = '))   # Escolhendo o spin;
n = int(input('n = '))

I = np.array([0, 0, 0])
II = np.array([a0, 0, 0])
III = np.array([0, a0, 0])
IV = np.array([a0, a0, 0])   # Definindo os vértices da celula unitária do crsital;
V = np.array([0, a0, a0])
VI = np.array([a0, a0, a0])
Z = np.array([a0/2, a0/2, a0/2])

rI = I - Z
rII = II - Z
rIII = III - Z   # Calculando a posição dos vértices em relação ao centro da célula;
rIV = IV - Z
rV = V - Z
rVI = VI - Z

q = np.linspace(-2*pi/a0, 2*pi/a0, 1000)   # Vetor de onda numérico;
Q = []

for i in q:
    Q.append(i)

dotI = []
dotII = []
dotIII = []
dotIV = []
dotV = []
dotVI = []

for j in Q:
    v = np.array([j, 0, 0])

    dotI.append(np.dot(v, rI))
    dotII.append(np.dot(v, rII))
    dotIII.append(np.dot(v, rIII))   # Produto interno numérico;
    dotIV.append(np.dot(v, rIV))
    dotV.append(np.dot(v, rV))
    dotVI.append(np.dot(v, rVI))

cos = []

for i in range(len(Q)):
    cos.append(np.cos(dotII[i]))

M = []

for i in range(len(Q)):
    M1 = np.zeros((n, n))

    for j in range(n):
        M1[j][j] = 4*J*S*(3-2*cos[i])

        if j == 0:
            M1[j][j + 1] = -4 * J * S * cos[i]

        if 1 <= j < n-1:
            M1[j][j + 1] = -4 * J * S * cos[i]
            M1[j][j - 1] = -4 * J * S * cos[i]

        if j == n-1:
            M1[j][j - 1] = -4 * J * S * cos[i]

    M.append(M1)

eig = []

for i in range(len(Q)):
    eigvals = np.linalg.eigvals(np.array(M[i]))

    eig.append(eigvals)

for i in range(n):
    E = []

    for j in range(len(Q)):
        e = eig[j][i]
        E.append(e)

    plt.scatter(q, E, color='r')

plt.show()
