from matplotlib.animation import FuncAnimation
import numpy as np
import matplotlib.pyplot as plt
'''В рамках борьбы с анархией сознания (ради порядка) начальные
параметры вводятся в СГС, а расчёты ведутся в (пс,мкм,N0)'''
Pi=3.141592653589
#ВВод начальных параметров, СГС
Nx=151#количество узлов сетки
Ny=31
Lx=5*10**-4# см 
Ly=1*10**-4# см 
N_0=1*10**11# см^-2
D_0=1*10**0# см^2/сек
R_a=1*10**-1# см^2/сек
V0D0_KbT=1*10**-10# см^4/сек
tau=500*10**-12# сек
#ДУ:z'=Az+Bz^2+C(лапласиан z)+D набла(z набла (z))
#z=n/N0, где n-концентрация
#НУ:n(x=0)=N0
#ГУ:n(y=0,Ly)=0
#Вычисление вспомогательных параметров, пересчёт в (пс,мкм,N0)
A=-1/tau*10**-12# пс^-1
B=-1*R_a*N_0*10**-12# пс^-1
C=D_0*10**-4#мкм^2/пс
D=V0D0_KbT*N_0*10**-4#мкм^2/пс
dx=Lx*10**4/(Nx-1)
dy=Ly*10**4/(Ny-1)#Шаги сетки по осям Х и У, мкм
t=0#время
dt=0.5#пс
Y, X = np.meshgrid(
    np.linspace(0, Ny-1, Ny),
    np.linspace(0, Nx-1, Nx)
)
Z=0*X*Y
#Элемент (x,y) определён как Z[x][y]
for j in range(0,Ny):
    Z[0][j]=1#НУ
    

def makeplot():
    cs = plt.contourf(dx*X,dy*Y,Z*N_0,levels=15)
    cbar=plt.colorbar(cs)
    cbar.set_label('Концентрация N, см^-2')
    plt.title('Время t = '+str(t)+' пс')
    plt.xlabel('Ось X, мкм')
    plt.ylabel('Ось Y, мкм')
def EulerStep(frame):
    global Z,t
    ZZ=0*X*Y
    
    for k in range(10):#количество итераций за фрейм
        for j in range(1,Ny-1):#Граничное условие
            N1=A*Z[0][j]+B*Z[0][j]*Z[0][j]
            N2=C*((Z[1][j]-2*Z[0][j])/(dx*dx)+(Z[0][j+1]-2*Z[0][j]+Z[0][j-1])/(dy*dy))
            N3=D*((Z[1][j]**2-2*Z[0][j]**2)/(2*dx*dx)+(Z[0][j+1]**2-2*Z[0][j]**2+Z[0][j-1]**2)/(2*dy*dy))
            ZZ[0][j]=Z[0][j]+(N1+N2+N3)*dt
        for i in range(1,Nx-1):#Граничное условие
            N1=A*Z[i][0]+B*Z[i][0]*Z[i][0]
            N2=C*((Z[i+1][0]-2*Z[i][0]+Z[i-1][0])/(dx*dx)+(Z[i][1]-2*Z[i][0])/(dy*dy))
            N3=D*((Z[i+1][0]**2-2*Z[i][0]**2+Z[i-1][0]**2)/(2*dx*dx)+(Z[i][1]**2-2*Z[i][0]**2)/(2*dy*dy))
            ZZ[i][0]=Z[i][0]+(N1+N2+N3)*dt
        for i in range(1,Nx-1):#Граничное условие
            N1=A*Z[i][Ny-1]+B*Z[i][Ny-1]*Z[i][Ny-1]
            N2=C*((Z[i+1][Ny-1]-2*Z[i][Ny-1]+Z[i-1][Ny-1])/(dx*dx)+(0-2*Z[i][Ny-1]+Z[i][Ny-2])/(dy*dy))
            N3=D*((Z[i+1][Ny-1]**2-2*Z[i][Ny-1]**2+Z[i-1][Ny-1]**2)/(2*dx*dx)+(0-2*Z[i][Ny-1]**2+Z[i][Ny-2]**2)/(2*dy*dy))
            ZZ[i][Ny-1]=Z[i][Ny-1]+(N1+N2+N3)*dt
        for i in range(1,Nx-1):
            for j in range(1,Ny-1):
                N1=A*Z[i][j]+B*Z[i][j]*Z[i][j]
                N2=C*((Z[i][j+1]-2*Z[i][j]+Z[i][j-1])/(dy*dy)+(Z[i+1][j]-2*Z[i][j]+Z[i-1][j])/(dx*dx))
                N3=D*((Z[i][j+1]**2-2*Z[i][j]**2+Z[i][j-1]**2)/(2*dy*dy)+(Z[i+1][j]**2-2*Z[i][j]**2+Z[i-1][j]**2)/(2*dx*dx))
                ZZ[i][j]=Z[i][j]+(N1+N2+N3)*dt
        t=t+dt
        Z=ZZ
    plt.clf()
    makeplot()
    if t==500:#Время (пс), на котором нужно остановить расчёт
        anim.event_source.stop()
plt.rcParams ['figure.figsize'] = [30*Lx/(Lx+Ly), 30*Ly/(Lx+Ly)]
fig,cs=plt.subplots()
makeplot()
#anim=FuncAnimation(fig,EulerStep,frames=None)
plt.show()
