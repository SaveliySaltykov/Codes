from matplotlib.animation import FuncAnimation
import numpy as np
import matplotlib.pyplot as plt
Pi=3.141592653589
Nx=41#количество узлов сетки
Ny=41
dx=2/(Nx-1)#сетка 2 на 2 мкм
dy=2/(Ny-1)
dt=1#пс
'''Система ДУ:
1)уравнение непрерывности N'=-div(N*vector(V))
2)ур-ие для Vx муторный крокодил
2)ур-ие для Vy'''
E=7.5
F=6.25
t=0#время
A=-2*10**-5#пс^-1
B=-5*10**-5#пс^-1
C=10**-6#мкм^2/пс
D=5*10**-6#мкм^2/пс
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

def D_x(N,i,j):#Взятие частной производной по x
    return (N[i+1][j]-N[i-1][j])/(2*dx)
def D_y(N,i,j):#Взятие частной производной по y
    return (N[i][j+1]-N[i][j-1])/(2*dy)
def Lap(N,i,j):#Взятие лапласиана
    x=(N[i+1][j]-2*N[i][j]+N[i-1][j])/(dx*dx)
    y=(N[i][j+1]-2*N[i][j]+N[i][j-1])/(dy*dy)
    return x+y

fig,cs=plt.subplots()
cs = plt.contourf(dx*(X-(Nx-1)/2),dy*(Y-(Ny-1)/2),Z*5*10**10,levels=15)
cbar=plt.colorbar(cs)
cbar.set_label('Концентрация N, см^-2')
plt.title('Время t = '+str(t)+' пс')
plt.xlabel('Ось X, мкм')
plt.ylabel('Ось Y, мкм')
def Eq1(Z,Vx,Vy,i,j):
    N1=-Z[i][j]*(D_x(Vx,i,j)+D_y(Vy,i,j))
    N2=-Vx[i][j]*D_x(Z,i,j)-Vy[i][j]*D_y(Z,i,j)
    return N1+N2
def Eq2(Z,Vx,Vy,i,j):
    N1=-Vx[i][j]*Eq1(Z,Vx,Vy,i,j)
    N2=-A*Z[i][j]*Vx[i][j]-B*D_x(Z,i,j)-C*D_x(Z,i,j)
    N3=D*Z[i][j]*Lap(Vx,i,j)
    N4=D*D_x(Z,i,j)*(D_x(Vx,i,j)-D_y(Vy,i,j))
    N5=D*D_y(Z,i,j)*(D_y(Vx,i,j)+D_x(Vy,i,j))
    return (N1+N2+N3+N4+N5)/Z[i,j]
def Eq3(Z,Vx,Vy,i,j):
    N1=-Vy[i][j]*Eq1(Z,Vx,Vy,i,j)
    N2=-A*Z[i][j]*Vy[i][j]-B*D_y(Z,i,j)-C*D_y(Z,i,j)
    N3=D*Z[i][j]*Lap(Vy,i,j)
    N4=D*D_y(Z,i,j)*(D_y(Vy,i,j)-D_x(Vx,i,j))
    N5=D*D_x(Z,i,j)*(D_x(Vy,i,j)+D_y(Vx,i,j))
    return (N1+N2+N3+N4+N5)/Z[i,j]
def EulerStep(frame):
    global Z,Vx,Vy,t
    ZZ=0*X*Y
    VVx=0*X*Y
    VVy=0*X*Y
    for k in range(10):#количество итераций за фрейм
        for i in range(1,Ny-1):#будет цикл от 1 до Nx-2
            for j in range(1,Nx-1):
                ZZ[i][j]=Z[i][j]+Eq1(Z,Vx,Vy,i,j)*dt#уравнение непрерывноести
                VVx[i][j]=Vx[i][j]+Eq2(Z,Vx,Vy,i,j)*dt#гидродинамическое
                VVy[i][j]=Vy[i][j]+Eq3(Z,Vx,Vy,i,j)*dt#гидродинамическое
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
    if t==200:#Время (пс), на котором нужно остановить расчёт
        anim.event_source.stop()
anim=FuncAnimation(fig,EulerStep,frames=None)
plt.show()
