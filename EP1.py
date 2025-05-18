#Autora: Carine Guzzi Santos 
#NUSP: 12556100
#Autora: Larissa 
#NUSP

""" 1. Resolva a equação, considerando um deslocamento inicial nulo 
e o sistema partindo do repouso. Use: """

# a) Método de Euler
# considerando o método de euler [y_{i+1}] = [y_i] + h[f(ti, [y_i])]
#imports
import matplotlib.pyplot as plt
import numpy as np


#variar o passo h utilizado em todo o código
h = 0.0025

# Constantes do sistema
M = 500
Ky = 20000
Cy = 700
tempo_total = 10  # tempo
n_passos = int(tempo_total / h)

# Condições iniciais
y1 = 0  # posição
y2 = 0  # velocidade
t = 0
yi = [y1, y2]

# Listas para armazenar resultados
tempos = []
posicoes = []
velocidades = []
aceleracoes = []

# Método de Euler
for i in range(n_passos):
    if t < 2:
        z = 0
    if t >= 2:
        z = 0.25

    y1 = yi[0]
    y2 = yi[1]
    y1pontoformula = y2
    y2pontoformula = (1/M)*((-Cy*(y2 - z))-(Ky*(y1 - z)))
    f = [y1pontoformula, y2pontoformula]

    yimaisum = [yi[0] + h*f[0], yi[1] + h*f[1]]
    yi = yimaisum
    t += h
    tempos.append(t)
    posicoes.append(yi[0])
    velocidades.append(yi[1])
    aceleracoes.append(y2pontoformula)

# Plotando os gráficos
plt.figure(figsize=(10, 8))

# Posição
plt.subplot(3, 1, 1)
plt.plot(tempos, posicoes, label="Posição y(t)", color='blue')
plt.title('Deslocamento')
plt.ylabel("Posição (m)")
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()

# Velocidade
plt.subplot(3, 1, 2)
plt.plot(tempos, velocidades, label="Velocidade y'(t)", color='green')
plt.title('Velocidade')
plt.ylabel("Velocidade (m/s)")
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()

# Aceleração
plt.subplot(3, 1, 3)
plt.plot(tempos, aceleracoes, label="Aceleração y''(t)", color='red')
plt.title('Aceleração')
plt.xlabel("Tempo (s)")
plt.ylabel("Aceleração (m/s²)")
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()

plt.tight_layout()
plt.show()


#b) método de Runge Kutta 4ªOrdem

# Parâmetros do sistema
t_max = 10     # tempo total de simulação

# Condições iniciais
y1 = 0        # posição inicial
y2 = 0        # velocidade inicial
t = 0

# Listas para armazenar resultados
T = [t]
Y1 = [y1]
Y2 = [y2]
A = [(1/M)*((-Cy*(y2 - 0)) - (Ky*(y1 - 0)))]  # aceleração inicial

# Função z(t) e z_ponto(t)
def z(t):
    return 0 if t < 2 else 0.25

def z_ponto(t):
    return 0  # degrau => derivada é zero exceto em t=2 (desconsiderado aqui)

# Função f(t, y1, y2)
def f(t, y1, y2):
    dz = z(t)
    dz_ponto = z_ponto(t)
    dy1 = y2
    dy2 = (1/M)*((-Cy*(y2 - dz_ponto)) - (Ky*(y1 - dz)))
    return dy1, dy2

# Método de Runge-Kutta de 4ª ordem
while t < t_max:
    k1_1, k1_2 = f(t, y1, y2)
    k2_1, k2_2 = f(t + h/2, y1 + h*k1_1/2, y2 + h*k1_2/2)
    k3_1, k3_2 = f(t + h/2, y1 + h*k2_1/2, y2 + h*k2_2/2)
    k4_1, k4_2 = f(t + h, y1 + h*k3_1, y2 + h*k3_2)

    y1 += h * (k1_1 + 2*k2_1 + 2*k3_1 + k4_1) / 6
    y2 += h * (k1_2 + 2*k2_2 + 2*k3_2 + k4_2) / 6
    t += h

    T.append(t)
    Y1.append(y1)
    Y2.append(y2)
    A.append((1/M)*((-Cy*(y2 - z_ponto(t))) - (Ky*(y1 - z(t)))))

# Plotar os resultados
plt.figure(figsize=(10, 12))

# Deslocamento
plt.subplot(3, 1, 1)
plt.plot(T, Y1, label='y(t) - Deslocamento', color='#1f77b4', linewidth=2)
plt.title('Deslocamento da Massa ao Longo do Tempo', fontsize=14)
plt.xlabel('Tempo (s)', fontsize=12)
plt.ylabel('Deslocamento (m)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=10)

# Velocidade
plt.subplot(3, 1, 2)
plt.plot(T, Y2, label='y\'(t) - Velocidade', color='#ff7f0e', linewidth=2)
plt.title('Velocidade da Massa ao Longo do Tempo', fontsize=14)
plt.xlabel('Tempo (s)', fontsize=12)
plt.ylabel('Velocidade (m/s)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=10)

# Aceleração
plt.subplot(3, 1, 3)
plt.plot(T, A, label='y\'\'(t) - Aceleração', color='#2ca02c', linewidth=2)
plt.title('Aceleração da Massa ao Longo do Tempo', fontsize=14)
plt.xlabel('Tempo (s)', fontsize=12)
plt.ylabel('Aceleração (m/s²)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=10)

# Ajusta layout
plt.tight_layout(pad=3.0)
plt.show()



#questão 2)
# Parâmetros do sistema
M = 500
Ky = 20000
Cy = 700
k = 10000
c = 350
L = 0.1

def z(t):
    return 0 if t < 2 else 0.25

def z_lombada(t, periodo=0.5):
    n = np.floor(2 / periodo)
    t1 = n * periodo
    t2 = (n + 10) * periodo
    if t < t1:
        return 0
    elif t1 <= t < t2:
        return 0.125 * (1 - np.cos(2 * np.pi * (t - t1) / periodo))
    else:
        return 0

def derivadas(t, y, y_ponto, u, u_ponto, m, tipo):
    if tipo == 'degrau':
        zt = z(t)
        zt_ponto = 0
    else:
        dt = 0.001
        zt = z_lombada(t)
        zt_ponto = (z_lombada(t + dt) - z_lombada(t - dt)) / (2 * dt)

    y2p = (-Ky * (y - zt) - Cy * (y_ponto - zt_ponto)) / M

    if m == 0:
        return y2p, 0

    u3 = u**3
    y2p += (k * u3 / (2 * L**2) + c * u_ponto) / (M + m)
    u2p = (-m * y2p - k * u3 / (2 * L**2) - c * u_ponto) / m
    return y2p, u2p

def runge_kutta(m, tempo, tipo):
    y = np.zeros(len(tempo))
    y_ponto = np.zeros(len(tempo))
    u = np.zeros(len(tempo))
    u_ponto = np.zeros(len(tempo))
    dt = tempo[1] - tempo[0]

    for i in range(len(tempo) - 1):
        k1_y, k1_u = derivadas(tempo[i], y[i], y_ponto[i], u[i], u_ponto[i], m, tipo)
        k2_y, k2_u = derivadas(tempo[i]+dt/2, y[i]+dt*y_ponto[i]/2, y_ponto[i]+dt*k1_y/2,
                               u[i]+dt*u_ponto[i]/2, u_ponto[i]+dt*k1_u/2, m, tipo)
        k3_y, k3_u = derivadas(tempo[i]+dt/2, y[i]+dt*y_ponto[i]/2, y_ponto[i]+dt*k2_y/2,
                               u[i]+dt*u_ponto[i]/2, u_ponto[i]+dt*k2_u/2, m, tipo)
        k4_y, k4_u = derivadas(tempo[i]+dt, y[i]+dt*y_ponto[i], y_ponto[i]+dt*k3_y,
                               u[i]+dt*u_ponto[i], u_ponto[i]+dt*k3_u, m, tipo)

        y_ponto[i+1] = y_ponto[i] + dt * (k1_y + 2*k2_y + 2*k3_y + k4_y) / 6
        y[i+1] = y[i] + dt * y_ponto[i]

        u_ponto[i+1] = u_ponto[i] + dt * (k1_u + 2*k2_u + 2*k3_u + k4_u) / 6
        u[i+1] = u[i] + dt * u_ponto[i]

    return y, u

#a)
# Gráficos y(t) e u(t) para degrau
tempo = np.arange(0, 5.0, 0.01)
epsilons = [0] + [e / 100 for e in np.arange(1, 20, 4.5)]

plt.figure(figsize=(12, 6))
for eps in epsilons:
    m = eps * M
    y, u = runge_kutta(m, tempo, tipo='degrau')
    plt.plot(tempo, y, label=f'y(t), ε = {eps:.3f}')
    if eps > 0:
        plt.plot(tempo, u, '--', label=f'u(t), ε = {eps:.3f}')
plt.xlabel('Tempo [s]')
plt.ylabel('Deslocamento [m]')
plt.title('Resposta ao degrau com e sem absorvedor')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()

#b)
# Gráficos y(t) para lombada
tempo = np.arange(0, 5.0, 0.001)
plt.figure(figsize=(12, 5))
for eps in epsilons:
    y, _ = runge_kutta(eps * M, tempo, tipo='lombada')
    plt.plot(tempo, y, label=f'ε = {eps * 100:.1f}%')
plt.title('Comparação de y(t) para diferentes valores de ε (lombada)')
plt.xlabel('Tempo [s]')
plt.ylabel('y(t) [m]')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()

# Gráficos u(t) para lombada
plt.figure(figsize=(12, 5))
for eps in epsilons[1:]:
    _, u = runge_kutta(eps * M, tempo, tipo='lombada')
    plt.plot(tempo, u, label=f'ε = {eps * 100:.1f}%')
plt.title('Comparação de u(t) para diferentes valores de ε (lombada)')
plt.xlabel('Tempo [s]')
plt.ylabel('u(t) [m]')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()

# Questão 3) Redução percentual da amplitude de y(t) em função de ε
faixa_eps = np.linspace(0.01, 0.19, 50)
max_y = []
for eps in faixa_eps:
    y, _ = runge_kutta(eps * M, tempo, tipo='lombada')
    max_y.append(np.max(np.abs(y)))

ref_y, _ = runge_kutta(0.0, tempo, tipo='lombada')
ref_max = np.max(np.abs(ref_y))
reducoes = [(1 - yi / ref_max) * 100 for yi in max_y]

plt.figure(figsize=(10, 5))
plt.plot(faixa_eps * 100, reducoes, 'o-', color='darkred')
plt.title('Redução percentual da amplitude de y(t) em função de ε')
plt.xlabel('ε [%]')
plt.ylabel('Redução de amplitude [%]')
plt.grid(True)
plt.tight_layout()
plt.show()


