from matplotlib.animation import FuncAnimation
import numpy as np
import matplotlib.pyplot as plt
Pi=3.141592653589
Nx=41#количество узлов сетки
Ny=41
dx=4/(Nx-1)#сетка 4 на 4 мкм
dy=4/(Ny-1)
dt=1#пс
'''Система ДУ:
1)уравнение непрерывности N'=-div(N*vector(V))
2)ур-ие для Vx муторный крокодил
2)ур-ие для Vy'''
tau_0=500#пс
tau_ex_ex=10**-1#пс
K_b=8.617333262*10**-5#эВ/K
T=50#K надо уточнять
m=0.067*0.511*0.5*10**9#эВ надо уточнять
N_0=500#мкм^-2
V_0=K_b*T*10**3#надо уточнять

E=7.5
F=6.25
t=0#время
A=1/tau_0#пс^-1
B=K_b*T/m#пс^-1
C=V_0*N_0/m#мкм^2/пс
D=tau_ex_ex*K_b*T#мкм^2/пс
#расстояние по сетке считается через мкм, время через пс
#z=n/N0, где N0=5*10^10 см^-2
Y, X = np.meshgrid(
    np.linspace(0, Ny-1, Ny),
    np.linspace(0, Nx-1, Nx)
)
Z=E/Pi*np.exp(-F*((dx*(X-(Nx-1)/2))**2+(dy*(Y-(Ny-1)/2))**2))
#Элемент (x,y) определён как Z[x][y]
Vx=0*X*Y#поле проекции скоростей Vx
Vy=0*X*Y#поле проекции скоростей Vy
'''Почему-то то, каким получится массив после умножения
зависит только от положения осей X и Y внутри meshgrid.
То есть Z=Y*X(=X*Y) даст такой же массив.'''

def D_x(N,x,y):#Взятие частной производной по x
    return (N[x+1][y]-N[x-1][y])/(2*dx)
def D_y(N,x,y):#Взятие частной производной по y
    return (N[x][y+1]-N[x][y-1])/(2*dy)
def Lap(N,x,y):#Взятие лапласиана
    N1=(N[x+1][y]-2*N[x][y]+N[x-1][y])/(dx*dx)
    N2=(N[x][y+1]-2*N[x][y]+N[x][y-1])/(dy*dy)
    return N1+N2

fig,cs=plt.subplots()
cs = plt.contourf(dx*(X-(Nx-1)/2),dy*(Y-(Ny-1)/2),Z*5*10**10,levels=15)
cbar=plt.colorbar(cs)
cbar.set_label('Концентрация N, см^-2')
plt.title('Время t = '+str(t)+' пс')
plt.xlabel('Ось X, мкм')
plt.ylabel('Ось Y, мкм')
def Eq1(x,y):
    N1=-Z[x][y]*(D_x(Vx,x,y)+D_y(Vy,x,y))
    N2=-Vx[x][y]*D_x(Z,x,y)-Vy[x][y]*D_y(Z,x,y)
    return N1+N2
def Eq2(x,y):
    N1=-Vx[x][y]*Eq1(x,y)
    N2=-A*Z[x][y]*Vx[x][y]-B*D_x(Z,x,y)-C*D_x(Z,x,y)
    N3=D*Z[x][y]*Lap(Vx,x,y)
    N4=D*D_x(Z,x,y)*(D_x(Vx,x,y)-D_y(Vy,x,y))
    N5=D*D_y(Z,x,y)*(D_y(Vx,x,y)+D_x(Vy,x,y))
    return (N1+N2+N3+N4+N5)/Z[x,y]
def Eq3(x,y):
    N1=-Vy[x][y]*Eq1(x,y)
    N2=-A*Z[x][y]*Vy[x][y]-B*D_y(Z,x,y)-C*D_y(Z,x,y)
    N3=D*Z[x][y]*Lap(Vy,x,y)
    N4=D*D_y(Z,x,y)*(D_y(Vy,x,y)-D_x(Vx,x,y))
    N5=D*D_x(Z,x,y)*(D_x(Vy,x,y)+D_y(Vx,x,y))
    return (N1+N2+N3+N4+N5)/Z[x,y]
def EulerStep(frame):
    global Z,Vx,Vy,t
    ZZ=0*X*Y
    VVx=0*X*Y
    VVy=0*X*Y
    for k in range(10):#количество итераций за фрейм анимации
        for i in range(1,Ny-1):#будет цикл от 1 до Nx-2
            for j in range(1,Nx-1):
                ZZ[i][j]=Z[i][j]+Eq1(i,j)*dt#уравнение непрерывноести
                VVx[i][j]=Vx[i][j]+Eq2(i,j)*dt#гидродинамическое
                VVy[i][j]=Vy[i][j]+Eq3(i,j)*dt#гидродинамическое
        t=t+dt
    plt.clf()
    Z=ZZ
    Vx=VVx
    Vy=VVy
    cs = plt.contourf(dx*(X-(Nx-1)/2),dy*(Y-(Ny-1)/2),Z*5*10**10,levels=15)
    cbar=plt.colorbar(cs)
    cbar.set_label('Концентрация N, см^-2')
    plt.title('Время t = '+str(t)+' пс')
    plt.xlabel('Ось X, мкм')
    plt.ylabel('Ось Y, мкм')
    if t==0:#Время (пс), на котором нужно остановить расчёт
        anim.event_source.stop()
#anim=FuncAnimation(fig,EulerStep,frames=None)
plt.show()
