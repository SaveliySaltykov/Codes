from matplotlib.animation import FuncAnimation
import numpy as np
import matplotlib.pyplot as plt
'''В рамках борьбы с анархией сознания (ради порядка) начальные
параметры вводятся в СГС, а расчёты ведутся в (пс,мкм,N0)'''
Pi=3.141592653589
#ВВод начальных параметров, СГС
Nx=41#количество узлов сетки
Ny=41
Lx=5*10**-4# см 
Ly=5*10**-4# см 
N_0=1*10**11# см^-2
D_0=3*10**0# см^2/сек
R_a=1*10**-1# см^2/сек
V0D0_KbT=1*10**-10# см^4/сек
tau=500*10**-12# сек
R_0=0.4*10**-4# см
#ДУ:z'=Az+Bz^2+C(лапласиан z)+D набла(z набла (z))
#z=n/N0, где n-концентрация
#НУ:z(r,0)=Exp(-(r/r0)^2)
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
Z=np.exp(-1/R_0/R_0*10**-8*((dx*(X-(Nx-1)/2))**2+(dy*(Y-(Ny-1)/2))**2))#+10**-10
#Элемент (x,y) определён как Z[x][y]
'''Почему-то то, каким получится массив после умножения
зависит только от положения осей X и Y внутри meshgrid.
То есть Z=Y*X(=X*Y) даст такой же массив.'''
def makeplot():
    cs = plt.contourf(dx*(X-(Nx-1)/2),dy*(Y-(Ny-1)/2),Z*N_0,levels=15)
    cbar=plt.colorbar(cs)
    cbar.set_label('Концентрация N, см^-2')
    plt.title('Время t = '+str(t)+' пс')
    plt.xlabel('Ось X, мкм')
    plt.ylabel('Ось Y, мкм')
def EulerStep(frame):
    global Z,t
    ZZ=0*X*Y
    for k in range(10):#количество итераций за фрейм
        for i in range(1,Nx-1):#будет цикл от 1 до Nx-2
            for j in range(1,Ny-1):
                N1=A*Z[i][j]+B*Z[i][j]*Z[i][j]
                N2=C*((Z[i][j+1]-2*Z[i][j]+Z[i][j-1])/(dy*dy)+(Z[i+1][j]-2*Z[i][j]+Z[i-1][j])/(dx*dx))
                N3=D*((Z[i][j+1]**2-2*Z[i][j]**2+Z[i][j-1]**2)/(2*dy*dy)+(Z[i+1][j]**2-2*Z[i][j]**2+Z[i-1][j]**2)/(2*dx*dx))
                ZZ[i][j]=Z[i][j]+(N1+N2+N3)*dt
        t=t+1
        Z=ZZ
    plt.clf()
    makeplot()   
    if t==1000:#Время (пс), на котором нужно остановить расчёт
        anim.event_source.stop()
fig,cs=plt.subplots()
makeplot()
anim=FuncAnimation(fig,EulerStep,frames=None)
plt.show()

    
    
