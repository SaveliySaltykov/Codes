from matplotlib.animation import FuncAnimation
import numpy as np
import matplotlib.pyplot as plt
Pi=3.141592653589
Nx=31#количество узлов сетки
Ny=31
dx=2/(Nx-1)#сетка 2 на 2 мкм
dy=2/(Ny-1)
#ДУ:z'=Az+Bz^2+C(лапласиан z)+D набла(z набла (z))
#расстояние по сетке считается через мкм, время через пс
#z=n/N0, где N0=5*10^10 см^-2
#НУ:z(r,0)=E/Pi*Exp(-Fr^2)
A=-2*10**-3#пс^-1
B=-5*10**-3#пс^-1
C=10**-4#мкм^2/пс
D=5*10**-4#мкм^2/пс
E=7.5
F=6.25
t=0#время

X, Y = np.meshgrid(
    np.linspace(0, Nx-1, Nx),
    np.linspace(0, Ny-1, Ny)
)
Z=E/Pi*np.exp(-F*((dx*(X-(Nx-1)/2))**2+(dy*(Y-(Ny-1)/2))**2))
#получилось так, что элемент (x,y) определён как Z[y][x]
fig,cs=plt.subplots()

cs = plt.contourf(dx*(X-(Nx-1)/2),dy*(Y-(Ny-1)/2),Z*5*10**10,levels=15)
cbar=plt.colorbar(cs)
cbar.set_label('Концентрация N, см^-2')
plt.title('Время t = '+str(t)+' пс')
plt.xlabel('Ось X, мкм')
plt.ylabel('Ось Y, мкм')

def EulerStep(frame):
    global Z,t
    ZZ=0*X*Y
    for k in range(10):#количество итераций за фрейм
        for i in range(1,Ny-1):#будет цикл от 1 до Nx-2
            for j in range(1,Nx-1):
                N1=A*Z[i][j]+B*Z[i][j]*Z[i][j]
                N2=C*((Z[i][j+1]-2*Z[i][j]+Z[i][j-1])/(dx*dx)+(Z[i+1][j]-2*Z[i][j]+Z[i-1][j])/(dy*dy))
                N3=D*((Z[i][j+1]**2-2*Z[i][j]**2+Z[i][j-1]**2)/(2*dx*dx)+(Z[i+1][j]**2-2*Z[i][j]**2+Z[i-1][j]**2)/(2*dy*dy))
                ZZ[i][j]=Z[i][j]+N1+N2+N3
            
        t=t+1
    plt.clf()
    Z=ZZ
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

    
    
