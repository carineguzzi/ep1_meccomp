#Autora: Carine Guzzi Santos 
#NUSP: 12556100
#Autora: Larissa 
#NUSP

""" 1. Resolva a equação, considerando um deslocamento inicial nulo 
e o sistema partindo do repouso. Use: """

#Equação M*dy2/dt2 + K(y-z) + C(dy/dt - dz/dt) = 0
# y1 = y
# y2 = yponto = y1ponto
# y1ponto = y2
# y2ponto = 

#Reescrevendo a equação: M*x1ponto + K(y-z) + C(dy/dt - dz/dt) = 0

# a) Método de Euler
# considerando o método de euler [y_{i+1}] = [y_i] + h[f(ti, [y_i])]
#constantes do sistema
import matplotlib.pyplot as plt

# Constantes do sistema
M = 500
Ky = 20000
Cy = 700
h = 0.0025
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
t_max = 5     # tempo total de simulação

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


