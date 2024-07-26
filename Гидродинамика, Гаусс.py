from matplotlib.animation import FuncAnimation
import numpy as np
import matplotlib.pyplot as plt
'''В рамках борьбы с анархией сознания (ради порядка) начальные
параметры вводятся в СГС, а расчёты ведутся в (пс,мкм,N0)'''
Pi=3.141592653589
#ВВод начальных параметров, СГС
Nx=100
Ny=100#Количество узлов сетки
Lx=8*10**-4# см 
Ly=8*10**-4# см 
tau_0=5.3*10**-13# сек 
KbT=5.52*10**-16# эрг
m=0.62*9.1*10**-28# г
V_0=1.21*1.6*10**-26# эрг*см^2
tau_ex_ex=5.6*10**-14# сек
N_0=1*10**11# см^-2
R_0=0.4*10**-4# см
'''z=n/N0, где n это 2D-концентрация
Система ДУ'ий:
1)уравнение непрерывности z'=-div(z*vector(V))
2)ур-ие для dVx/dt это муторный крокодил вида:
Vx'=...-A(...)-B(...)-C(...)+D(...)
3)ур-ие для Vy'''
#Вычисление вспомогательных параметров, пересчёт в (пс,мкм,N0)
A=1/tau_0*10**-12# пс^-1, 1сек=10^12пс
B=KbT/m*10**-16#(мкм/пс)^2, 1см=10^4мкм
C=V_0*N_0/m*10**-16# (мкм/пс)^2
D=tau_ex_ex*KbT/m*10**-4# мкм^2/пс
dx=Lx*10**4/(Nx-1)
dy=Ly*10**4/(Ny-1)#Шаги сетки по осям Х и У, мкм
t=0#время, пс
dt=1#пс
Y, X = np.meshgrid(
    np.linspace(0, Ny-1, Ny),
    np.linspace(0, Nx-1, Nx)
)
Z=np.exp(-1/R_0/R_0*10**-8*((dx*(X-(Nx-1)/2))**2+(dy*(Y-(Ny-1)/2))**2))
for i in range(0,Nx):
    for j in range(0,Ny):
        if Z[i][j]<10**-20:
            Z[i][j]=0
#Элемент (x,y) определён как Z[x][y]
Vx=0*X*Y#+dx*(X-(Nx-1)/2)/R_0*np.sqrt(KbT/m)*10**-12# мкм/пс
#поле проекции скоростей Vx
Vy=0*X*Y#+dy*(Y-(Ny-1)/2)/R_0*np.sqrt(KbT/m)*10**-12# мкм/пс
#поле проекции скоростей Vy
'''То, каким получится массив после умножения
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
def makeplot():
    cs = plt.contourf(dx*(X-(Nx-1)/2),dy*(Y-(Ny-1)/2),Z*N_0,levels=15)
    cbar=plt.colorbar(cs)
    cbar.set_label('Концентрация N, см^-2')
    plt.title('Время t = '+str(t)+' пс (Гидродинамика)')
    plt.xlabel('Ось X, мкм')
    plt.ylabel('Ось Y, мкм')
def Eq1(x,y):
    N1=-Z[x][y]*(D_x(Vx,x,y)+D_y(Vy,x,y))
    N2=-Vx[x][y]*D_x(Z,x,y)-Vy[x][y]*D_y(Z,x,y)
    return N1+N2
def Eq2(x,y):
    N1=-Vx[x][y]*Eq1(x,y)/Z[x,y]
    N2=-A*Z[x][y]*Vx[x][y]-B*D_x(Z,x,y)/Z[x,y]-C*D_x(Z,x,y)
    N3=D*Z[x][y]*Lap(Vx,x,y)
    N4=D*D_x(Z,x,y)*(D_x(Vx,x,y)-D_y(Vy,x,y))/Z[x,y]
    N5=D*D_y(Z,x,y)*(D_y(Vx,x,y)+D_x(Vy,x,y))/Z[x,y]
    return (N1+N2+N3+N4+N5)
def Eq3(x,y):
    N1=-Vy[x][y]*Eq1(x,y)/Z[x,y]
    N2=-A*Z[x][y]*Vy[x][y]-B*D_y(Z,x,y)/Z[x,y]-C*D_y(Z,x,y)
    N3=D*Z[x][y]*Lap(Vy,x,y)
    N4=D*D_y(Z,x,y)*(D_y(Vy,x,y)-D_x(Vx,x,y))/Z[x,y]
    N5=D*D_x(Z,x,y)*(D_x(Vy,x,y)+D_y(Vx,x,y))/Z[x,y]
    return (N1+N2+N3+N4+N5)
def EulerStep(frame):
    global Z,Vx,Vy,t
    plt.clf()
    makeplot()
    if round(t)==1000:#Время (пс), на котором нужно остановить расчёт
        anim.event_source.stop()
    ZZ=0*X*Y
    VVx=0*X*Y
    VVy=0*X*Y
    for k in range(10):#количество итераций за фрейм анимации
        for i in range(1,Nx-1):#цикл от 1 до Nx-2
            for j in range(1,Ny-1):
                ZZ[i][j]=Z[i][j]+Eq1(i,j)*dt#уравнение непрерывноести
                if Z[i][j]==0:
                    VVx[i][j]=Vx[i][j]
                    VVy[i][j]=Vy[i][j]
                else:
                    VVx[i][j]=Vx[i][j]+Eq2(i,j)*dt#гидродинамическое
                    VVy[i][j]=Vy[i][j]+Eq3(i,j)*dt#гидродинамическое
        """for i in range(1,Nx-1):
            VVy[i][0]=VVy[i][1]
            VVy[i][Ny-1]=VVy[i][Ny-2]
            VVx[i][0]=VVx[i][1]
            VVx[i][Ny-1]=VVx[i][Ny-2]
        for j in range(1,Nx-1):
            VVx[0][j]=VVx[1][j]
            VVx[Nx-1][j]=VVx[Nx-2][j]
            VVy[0][j]=VVy[1][j]
            VVy[Nx-1][j]=VVy[Nx-2][j]"""
        Z=ZZ
        Vx=VVx
        Vy=VVy
        t=round(t,2)+dt
fig,cs=plt.subplots()
makeplot()
anim=FuncAnimation(fig,EulerStep,frames=None)

plt.show()
print((np.sqrt(KbT/m))*10**-8)


