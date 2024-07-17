from matplotlib.animation import FuncAnimation
import numpy as np
import matplotlib.pyplot as plt
Pi=3.141592653589
Nx=121#количество узлов сетки
Ny=31
dx=8/(Nx-1)#сетка 2 на 8 мкм
dy=2/(Ny-1)
plt.rcParams ['figure.figsize'] = [16, 4]
#ДУ:z'=Az+Bz^2+C(лапласиан z)+D набла(z набла (z))
#расстояние по сетке считается через мкм, время через пс
#z=n/N0, где N0=5*10^10 см^-2
#НУ:z(x=0)=E/Pi
#ГУ:на стенках z=0 набла (z)=0
A=-2*10**-3#пс^-1
B=-5*10**-3#пс^-1
C=10**-4#мкм^2/пс
D=5*10**-4#мкм^2/пс
E=7.5
t=0#время, пс
X, Y = np.meshgrid(
    np.linspace(0, Nx-1, Nx),
    np.linspace(0, Ny-1, Ny)
)
Z=0*X*Y
#получилось так, что элемент (x,y) определён как Z[y][x]
for i in range(1,Ny-1):
    Z[i][0]=E/Pi#НУ
    
fig,cs=plt.subplots()
cs = plt.contourf(dx*X,dy*Y,Z*5*10**10,levels=15)
cbar=plt.colorbar(cs)
cbar.set_label('Концентрация N, см^-2')
plt.title('Время t = '+str(t)+' пс')
plt.xlabel('Ось X, мкм')
plt.ylabel('Ось Y, мкм')
def EulerStep(frame):
    global Z,t
    ZZ=0*X*Y
    
    for k in range(1):#количество итераций за фрейм
        for i in range(1,Ny-1):#Граничное условие
            N1=A*Z[i][0]+B*Z[i][0]*Z[i][0]
            N2=C*((Z[i][1]-2*Z[i][0])/(dx*dx)+(Z[i+1][0]-2*Z[i][0]+Z[i-1][0])/(dy*dy))
            N3=D*((Z[i][1]**2-2*Z[i][0]**2)/(2*dx*dx)+(Z[i+1][0]**2-2*Z[i][0]**2+Z[i-1][0]**2)/(2*dy*dy))
            ZZ[i][0]=Z[i][0]+N1+N2+N3
        for j in range(1,Nx-1):#Граничное условие
            N1=A*Z[0][j]+B*Z[0][j]*Z[0][j]
            N2=C*((Z[0][j+1]-2*Z[0][j]+Z[0][j-1])/(dx*dx)+(Z[1][j]-2*Z[0][j])/(dy*dy))
            N3=D*((Z[0][j+1]**2-2*Z[0][j]**2+Z[0][j-1]**2)/(2*dx*dx)+(Z[1][j]**2-2*Z[0][j]**2)/(2*dy*dy))
            ZZ[0][j]=Z[0][j]+N1+N2+N3
        for j in range(1,Nx-1):#Граничное условие
            N1=A*Z[Ny-1][j]+B*Z[Ny-1][j]*Z[Ny-1][j]
            N2=C*((Z[Ny-1][j+1]-2*Z[Ny-1][j]+Z[Ny-1][j-1])/(dx*dx)+(0-2*Z[Ny-1][j]+Z[Ny-2][j])/(dy*dy))
            N3=D*((Z[Ny-1][j+1]**2-2*Z[Ny-1][j]**2+Z[Ny-1][j-1]**2)/(2*dx*dx)+(0-2*Z[Ny-1][j]**2+Z[Ny-2][j]**2)/(2*dy*dy))
            ZZ[Ny-1][j]=Z[Ny-1][j]+N1+N2+N3
        for i in range(1,Ny-1):
            for j in range(1,Nx-1):
                N1=A*Z[i][j]+B*Z[i][j]*Z[i][j]
                N2=C*((Z[i][j+1]-2*Z[i][j]+Z[i][j-1])/(dx*dx)+(Z[i+1][j]-2*Z[i][j]+Z[i-1][j])/(dy*dy))
                N3=D*((Z[i][j+1]**2-2*Z[i][j]**2+Z[i][j-1]**2)/(2*dx*dx)+(Z[i+1][j]**2-2*Z[i][j]**2+Z[i-1][j]**2)/(2*dy*dy))
                ZZ[i][j]=Z[i][j]+N1+N2+N3
        t=t+1
        
    plt.clf()
    Z=ZZ
    cs = plt.contourf(dx*X,dy*Y,Z*5*10**10,levels=15)
    cbar=plt.colorbar(cs)
    cbar.set_label('Концентрация N, см^-2')
    plt.title('Время t = '+str(t)+' пс')
    plt.xlabel('Ось X, мкм')
    plt.ylabel('Ось Y, мкм')
anim=FuncAnimation(fig,EulerStep,frames=None)
plt.show()
