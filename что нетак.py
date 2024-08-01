from matplotlib.animation import FuncAnimation
import numpy as np
import matplotlib.pyplot as plt

'''В рамках борьбы с анархией сознания (ради порядка) начальные
параметры вводятся в СГС, а расчёты ведутся в (пс,мкм,N0)'''
Pi=3.141592653589
#ВВод начальных параметров, СГС
Nx=33#количество узлов сетки
Ny=19
Lx=3*10**-4# см 
Ly=0.5*10**-4# см 
tau_0=5.3*10**-12# сек 
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
#Вычисление вспомогательных параметров, пересчёт в (пс,мкм,?N0?возможно, надо отказаться от замены z=n/N0)
A=1/tau_0*10**-12# пс^-1, 1сек=10^12пс
B=KbT/m*10**-16#(мкм/пс)^2, 1см=10^4мкм
C=V_0*N_0/m*10**-16# (мкм/пс)^2
D=tau_ex_ex*KbT/m*10**-4# мкм^2/пс
dx=Lx*10**4/(Nx-3)
dy=Ly*10**4/(Ny-3)#Шаги сетки по осям Х и У, мкм
t=0#время, пс
dt=10**-10#пс
Y, X = np.meshgrid(
    np.linspace(0, Ny-1, Ny),
    np.linspace(0, Nx-1, Nx)
)
Z=0*X*Y+10**-10
#Элемент (x,y) определён как Z[x][y]
for j in range(2,Ny-2):
    Z[1][j]=1
Vx=0*X*Y#+dx*(X-(Nx-1)/2)/R_0*np.sqrt(KbT/m)*10**-12# мкм/пс
#поле проекции скоростей Vx
Vy=0*X*Y#+dy*(Y-(Ny-1)/2)/R_0*np.sqrt(KbT/m)*10**-12# мкм/пс
#поле проекции скоростей Vy
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
def makeplot():
    XLine=np.linspace(0, Nx-3, Nx-2)
    Yline=np.linspace(0, Ny-3, Ny-2)
    YY, XX = np.meshgrid(
        Yline,
        XLine
    )
    ZZ=0*YY*XX
    VVx=0*YY*XX
    VVy=0*YY*XX
    for i in range(0,Nx-2):
        for j in range(0,Ny-2):
            ZZ[i][j]=Z[i+1][j+1]
    for i in range(0,Nx-2):
        for j in range(0,Ny-2):
            VVx[i][j]=Vx[i+1][j+1]
    for i in range(0,Nx-2):
        for j in range(0,Ny-2):
            VVy[i][j]=Vy[i+1][j+1]
    #Пробую интерполяцию для сглаживания
    """YYY, XXX = np.meshgrid(
        np.linspace(0, (Ny-2)*2-2, (Ny-2)*2-1),
        np.linspace(0, (Nx-2)*2-2, (Nx-2)*2-1)
    )
    points=[[2*x,2*y] for x in XLine for y in Yline]

    for i in range(0,Nx-2):
        ZZ[i][0]=ZZ[i][1]
        ZZ[i][Nx-3]=ZZ[i][Nx-4]
    values=np.append([],ZZ)
    
    ZZZ = griddata(points,values,(XXX,YYY),method='cubic')
    for i in range(0,2*Nx-5):
        for j in range (0,2*Ny-5):
            if ZZZ[i][j]<0:
                ZZZ[i][j]=0
    for i in range(0,2*Nx-5):
        for j in range(0,Ny-2):
            ZZZ[i][j]=ZZZ[i][(2*Ny-6)-j]"""
    cs = plt.contourf(dx*XX,dy*YY,ZZ*N_0,cmap='jet',levels=15)#ZZ*N_0
    cbar=plt.colorbar(cs)
    cbar.set_label('Концентрация N, см^-2') 
    plt.title('Время t = '+str(t)+' пс (Гидродинамика)')
    if t>100:
        plt.title('Время t = '+str(round(t))+' пс (Гидродинамика)')
    plt.xlabel('Ось X, мкм')
    plt.ylabel('Ось Y, мкм')
def Eq1(N,NVx,NVy,x,y):
    N1=-N[x][y]*(D_x(NVx,x,y)+D_y(NVy,x,y))
    N2=-NVx[x][y]*D_x(N,x,y)-NVy[x][y]*D_y(N,x,y)
    return N1+N2
def Eq2(N,NVx,NVy,x,y):
    N=N+10**-30
    N1=NVx[x][y]*((D_x(NVx,x,y)+D_y(NVy,x,y))+NVx[x][y]*D_x(np.log(N),x,y)+NVy[x][y]*D_y(np.log(N),x,y))
    N2=-A*NVx[x][y]-B*D_x(np.log(N),x,y)-C*D_x(N,x,y)
    N3=D*Lap(NVx,x,y)#*Z[x][y]
    N4=D*D_x(np.log(N),x,y)*(D_x(NVx,x,y)-D_y(NVy,x,y))
    N5=D*D_y(np.log(N),x,y)*(D_y(NVx,x,y)+D_x(NVy,x,y))
    return N1+N2+N3+N4+N5
def Eq3(N,NVx,NVy,x,y):
    N=N+10**-30
    N1=NVy[x][y]*((D_x(NVx,x,y)+D_y(NVy,x,y))+NVx[x][y]*D_x(np.log(N),x,y)+NVy[x][y]*D_y(np.log(N),x,y))
    N2=-A*NVy[x][y]-B*D_y(np.log(N),x,y)-C*D_y(N,x,y)
    N3=D*Lap(NVy,x,y)#*Z[x][y]
    N4=D*D_y(np.log(N),x,y)*(D_y(NVy,x,y)-D_x(NVx,x,y))
    N5=D*D_x(np.log(N),x,y)*(D_x(NVy,x,y)+D_y(NVx,x,y))
    return N1+N2+N3+N4+N5

   
    
def Granich(N,NVx,NVy):
    for i in range(0,Nx):
        N[i][0]=N[i][1]
        N[i][Ny-1]=N[i][Ny-2]
        NVy[i][0]=-NVy[i][1]
        NVy[i][Ny-1]=-NVy[i][Ny-2]
        NVx[i][0]==-NVx[i][1]
        NVx[i][Ny-1]=-NVx[i][Ny-2]
    for j in range(0,Ny):
        N[0][j]=N[1][j]
        NVx[0][j]=-NVx[1][j]
        NVy[0][j]=-NVy[1][j]
        N[Nx-1][j]=N[Nx-2][j]
        NVx[Nx-1][j]=-NVx[Nx-2][j]
        NVy[Nx-1][j]=-NVy[Nx-2][j]
    #N[1][1]=N[2][1]
    #N[Nx-2][1]=N[Nx-3][1]
    #N[1][Ny-2]=N[2][Ny-2]
    #N[Nx-2][Ny-2]=N[Nx-3][Ny-2]
    return N,NVx,NVy    
P=1
def EulerStep(frame):
    global Z,Vx,Vy,t,dt,P
    plt.clf()
    makeplot()
    
    DZ1=0*X*Y
    DVx1=0*X*Y
    DVy1=0*X*Y
    DZ2=0*X*Y
    DVx2=0*X*Y
    DVy2=0*X*Y
    DZ3=0*X*Y
    DVx3=0*X*Y
    DVy3=0*X*Y
    DZ4=0*X*Y
    DVx4=0*X*Y
    DVy4=0*X*Y
    
    for k in range(P):#количество итераций за фрейм анимации
        for i in range(1,Nx-1):
            Vx[i][1]=0
            Vx[i][Ny-2]=0
            Vy[i][1]=0
            Vy[i][Ny-2]=0
        for j in range(1,Ny-1):
            Vy[1][j]=0
            Vy[Nx-2][j]=0
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
            dt=0.25
            P=10
        if t>100:
            P=100
        
        """if t>10*dt:
            dt=dt*10
        if t>0.01:
            dt=0.001
        if t>0.2:
            dt=0.01
        if t>2:
            dt=0.1
        if round(t)==250:
            plt.clf()
            makeplot()
        if round(t)==500:
            plt.clf()
            makeplot()  """  
        if round(t)==1000:#Время (пс), на котором нужно остановить расчёт
            plt.clf()
            makeplot()
            plt.title('Время t = '+str(round(t))+' пс')
            anim.event_source.stop()
        
        """for i in range(0,Nx):#Почему-то эти условия для слоёв за стенкой вызывают рост частиц полсе первой пикосекунды
            Z[i][0]=2*Z[i][1]-Z[i][2]
            Z[i][Ny-1]=2*Z[i][Ny-2]-Z[i][Ny-3]
            Vy[i][0]=2*Vy[i][1]-Vy[i][2]
            Vy[i][Ny-1]=2*Vy[i][Ny-2]-Vy[i][Ny-3]
            Vx[i][0]==2*Vx[i][1]-Vx[i][2]
            Vx[i][Ny-1]=2*Vx[i][Ny-2]-Vx[i][Ny-3]
        for j in range(0,Ny):
            Z[0][j]=2*Z[1][j]-Z[2][j]
            Vx[0][j]=2*Vx[1][j]-Vx[2][j]
            Vy[0][j]=2*Vy[1][j]-Vy[2][j]
            Z[Nx-1][j]=2*Z[Nx-2][j]-Z[Nx-3][j]
            Vx[Nx-1][j]=2*Vx[Nx-2][j]-Vx[Nx-3][j]
            Vy[Nx-1][j]=2*Vy[Nx-2][j]-Vy[Nx-3][j]"""
        Z,Vx,Vy=Granich(Z,Vx,Vy)
        for i in range(1,Nx-1): 
            for j in range(1,Ny-1):
                DZ1[i][j]=Eq1(Z,Vx,Vy,i,j)#уравнение непрерывноести
                DVx1[i][j]=Eq2(Z,Vx,Vy,i,j)#гидродинамическое
              #  DVy1[i][j]=Eq3(Z,Vx,Vy,i,j)#гидродинамическое
        DZ1,DVx1,DVy1=Granich(DZ1,DVx1,DVy1)
        for i in range(1,Nx-1): 
            for j in range(1,Ny-1):
                DZ2[i][j]=Eq1(Z+dt*DZ1/2,Vx+dt*DVx1/2,Vy+dt*DVy1/2,i,j)#уравнение непрерывноести
                DVx2[i][j]=Eq2(Z+dt*DZ1/2,Vx+dt*DVx1/2,Vy+dt*DVy1/2,i,j)#гидродинамическое
             #   DVy2[i][j]=Eq3(Z+dt*DZ1/2,Vx+dt*DVx1/2,Vy+dt*DVy1/2,i,j)#гидродинамическое
        DZ2,DVx2,DVy2=Granich(DZ2,DVx2,DVy2)
        for i in range(1,Nx-1): 
            for j in range(1,Ny-1):
                DZ3[i][j]=Eq1(Z+dt*DZ2/2,Vx+dt*DVx2/2,Vy+dt*DVy2/2,i,j)#уравнение непрерывноести
                DVx3[i][j]=Eq2(Z+dt*DZ2/2,Vx+dt*DVx2/2,Vy+dt*DVy2/2,i,j)#гидродинамическое
              #  DVy3[i][j]=Eq3(Z+dt*DZ2/2,Vx+dt*DVx2/2,Vy+dt*DVy2/2,i,j)#гидродинамическое
        DZ3,DVx3,DVy3=Granich(DZ3,DVx3,DVy3)
        for i in range(1,Nx-1): 
            for j in range(1,Ny-1):
                DZ1[i][j]=Eq1(Z+dt*DZ3,Vx+dt*DVx3,Vy+dt*DVy3,i,j)#уравнение непрерывноести
                DVx1[i][j]=Eq2(Z+dt*DZ3,Vx+dt*DVx3,Vy+dt*DVy3,i,j)#гидродинамическое
              #  DVy1[i][j]=Eq3(Z+dt*DZ3,Vx+dt*DVx3,Vy+dt*DVy3,i,j)#гидродинамическое
        DZ4,DVx4,DVy4=Granich(DZ4,DVx4,DVy4)
        Z=Z+dt*(DZ1+2*DZ2+2*DZ3+DZ4)/6
        Vx=Vx+dt*(DVx1+2*DVx2+2*DVx3+DVx4)/6
        Vy=0*X*Y#Vy+dt*(DVy1+2*DVy2+2*DVy3+DVy4)/6
        Count=0
        for i in range(1,Nx-1): 
            for j in range(1,Ny-1):
                Count=Count+Z[i][j]
        if t<100:
            print(str(Count)+'\t'+str(Z[2][2]))
        t=t+dt
        
plt.rcParams ['figure.figsize'] = [30*Lx/(Lx+Ly), 30*Ly/(Lx+Ly)]
fig,cs=plt.subplots()
makeplot()
anim=FuncAnimation(fig,EulerStep,frames=None)
plt.show()
