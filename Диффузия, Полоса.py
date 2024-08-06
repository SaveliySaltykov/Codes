from matplotlib.animation import FuncAnimation
import numpy as np
import matplotlib.pyplot as plt
'''В рамках борьбы с анархией сознания (ради порядка) начальные
параметры вводятся в СГС, а расчёты ведутся в (пс,мкм,N0)'''
Pi=3.141592653589
#ВВод начальных параметров, СГС
Nx=33#количество узлов сетки
Ny=20
Lx=3*10**-4# см 
Ly=1.0*10**-4# см 
N_0=1*10**11# см^-2
D_0=1*10**0# см^2/сек
R_a=1*10**-1# см^2/сек
V0D0_KbT=1*10**-10# см^4/сек
tau=500*10**-12# сек
R_0=0.4*10**-4# см
#ДУ:z'=Az+Bz^2+C(лапласиан z)+D набла(z набла (z))
#z=n/N0, где n-концентрация
#НУ:n(x=0)=N0
#ГУ:n(y=0,Ly)=0
#Вычисление вспомогательных параметров, пересчёт в (пс,мкм,N0)
A=-1/tau*10**-12# пс^-1
B=-1*R_a*N_0*10**-12# пс^-1
C=D_0*10**-4#мкм^2/пс
D=V0D0_KbT*N_0*10**-4#мкм^2/пс
dx=Lx*10**4/(Nx-3)
dy=Ly*10**4/(Ny-3)#Шаги сетки по осям Х и У, мкм
t=0#время
dt=10**-3#пс
Y, X = np.meshgrid(
    np.linspace(0, Ny-1, Ny),
    np.linspace(0, Nx-1, Nx)
)
Z=np.exp(-1/R_0/R_0*10**-8*((dx*(X-1))**2+(dy*(Y-(Ny-1)/2))**2))+10**-10
#Элемент (x,y) определён как Z[x][y]

    

def makeplot():
    YY, XX = np.meshgrid(
        np.linspace(0, Ny-3, Ny-2),
        np.linspace(0, Nx-3, Nx-2)
    )
    ZZ=0*YY*XX
    for i in range(0,Nx-2):
        for j in range(0,Ny-2):
            if t<300:
                ZZ[i][j]=Z[i+1][j+1]
            else:
                ZZ[i][j]=Z[i+1][Ny//2]
    
    
    cs = plt.contourf(dx*XX,dy*YY,ZZ*N_0,cmap='jet',levels=15)#ZZ*N_0
    cbar=plt.colorbar(cs)
    cbar.set_label('Концентрация N, см^-2') 
    plt.title('Время t = '+str(t)+' пс (Диффузия)')
    if t>100:
        plt.title('Время t = '+str(round(t))+' пс (Диффузия)')
    plt.xlabel('Ось X, мкм')
    plt.ylabel('Ось Y, мкм')
P=10
def EulerStep(frame):
    global Z,t,dt,P
    plt.clf()
    makeplot()
    ZZ=0*X*Y
    
    for k in range(P):#количество итераций за фрейм
        if round(t)==250:
            plt.clf()
            makeplot()
            plt.savefig('DP1.png')
        if round(t)==500:
            plt.clf()
            makeplot()
            plt.savefig('DP2.png')
        if round(t)==750:
            plt.clf()
            makeplot()
            plt.savefig('DP3.png')
        if round(t)==1000:
            plt.clf()
            makeplot()
            plt.savefig('DP4.png')
            anim.event_source.stop()
        if t>10*dt:
            dt=dt*10
        if t>0.001:
            dt=0.0001
        if t>0.01:
            dt=0.001
        if t>0.1:
            dt=0.01
        if t>1:
            dt=0.1
        if t>10:
            dt=0.5
            P=10
        if t>100:
            P=100
        for j in range(1,Ny-1):#Граничное условие
            Z[0][j]=Z[1][j]
            Z[Nx-1][j]=Z[Nx-2][j]
        for i in range(1,Nx-1):#Граничное условие
            Z[i][0]=Z[i][1]
            Z[i][Ny-1]=Z[i][Ny-2]
        for i in range(1,Nx-1):
            for j in range(1,Ny-1):
                N1=A*Z[i][j]+B*Z[i][j]*Z[i][j]
                N2=C*((Z[i][j+1]-2*Z[i][j]+Z[i][j-1])/(dy*dy)+(Z[i+1][j]-2*Z[i][j]+Z[i-1][j])/(dx*dx))
                N3=D*((Z[i][j+1]**2-2*Z[i][j]**2+Z[i][j-1]**2)/(2*dy*dy)+(Z[i+1][j]**2-2*Z[i][j]**2+Z[i-1][j]**2)/(2*dx*dx))
                ZZ[i][j]=Z[i][j]+(N1+N2+N3)*dt 
        t=t+dt
        Z=ZZ
    
    
plt.rcParams ['figure.figsize'] = [27*Lx/(Lx+Ly), 27*Ly/(Lx+Ly)]
fig,cs=plt.subplots()
makeplot()
plt.savefig('DP0.png')
anim=FuncAnimation(fig,EulerStep,frames=None)
plt.show()
