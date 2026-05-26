import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import matplotlib.animation as animation

fig=plt.figure()
fig.set_dpi(100)
fig.set_size_inches(7,6.5)

#Datos a introducir:
L1=0.1      #longitud del péndulo 1
L2=0.1      #longitud del pendulo 2
m1=1.0      #masa del pendulo 1
m2=1.0      #masa del pendulo 2
g=9.81      #aceleracion de la gravedad
tf=10.0     #tiempo de simulacion

m12=m1+m2
par=[L1,L2,m1,m2,m12,g]

#Ecuaciones del movimiento del pendulo doble:
def double_pendulum(z,t,par):
    z1,z2,z3,z4=z  

    sinz=np.sin(z1-z2)
    cosz=np.cos(z1-z2)
    sinz1=np.sin(z1)
    sinz2=np.sin(z2)
    z42=z4*z4
    z32=z3*z3
    coszsq=cosz*cosz

    dzdt=[z3,z4,(-m2*L1*z32*sinz*cosz+g*m2*sinz2*cosz-m2*L2*z42*sinz-m12*g*sinz1)/(L1*m12-m2*L1*coszsq),
         (m2*L2*z42*sinz*cosz+g*sinz1*cosz*m12+L1*z32*sinz*m12-g*sinz2*m12)/(L2*m12-m2*L2*coszsq)]
    return dzdt


#Llamada a odeint que resuelve las ecuaciones de movimiento:
nt=500  #Numero de intervalos de tiempo
dt=tf/nt

#Angulos iniciales con respecto a la vertical de cada pendulo:
#Angulos iniciales en grados:
theta1_0_deg=20.0
theta2_0_deg=20.0
#Angulos iniciales en radianes:
theta1_0=theta1_0_deg*np.pi/180.0
theta2_0=theta2_0_deg*np.pi/180.0
z0=[theta1_0,theta2_0,0.0,0.0] #Valores iniciales

#z0 = [θ10, θ20, ω10, ω20] #Consideremos que no tienen ω0
t=np.linspace(0,tf,nt)
abserr = 1.0e-8
relerr = 1.0e-6
z=odeint(double_pendulum,z0,t,args=(par,),atol=abserr, rtol=relerr)
plt.close('all')
"""
#Angulos en funcion del tiempo:
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.set_xlabel('tiempo (s)')
ax1.set_ylabel('$\\theta_1$(rad)', color='b', fontsize=15)
ax2.set_ylabel('$\\theta_2$(rad)', color='r', fontsize=15)
plt.title('Ángulos en función del tiempo')
ax1.tick_params('y', colors='b')
ax2.tick_params('y', colors='r')
ax1.set_xlim(xmin=0.,xmax=tf) #limites del eje x
line1, = ax1.plot(t[:],z[:,0], linewidth=2, color='b')
line2, = ax1.plot(t[:],z[:,1], linewidth=2, color='r')
plt.show()

#Espacio de fase ω(θ):
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.set_xlabel('$\\theta$(rad)')
ax1.set_ylabel('$\dot{\\theta}_1$(rad/s)', color='b', fontsize=15)
ax2.set_ylabel('$\dot{\\theta}_2$(rad/s)', color='r', fontsize=15)
plt.title('Espacio de Fase')
ax1.tick_params('y', colors='b')
ax2.tick_params('y', colors='r')
line1, = plt.plot(z[:,0], z[:,2], linewidth=2, c='b')
line2, = plt.plot(z[:,1], z[:,3], linewidth=2, c='r')
plt.show()

#Energia de cada masa y del sistema:
E1 = 0.5*m1*z[:,2]**2*L1**2 - m1*g*L1*np.cos(z[:,0])
E2 = 0.5*m2*(z[:,2]**2*L1**2 + z[:,3]**2*L2**2 + 2*z[:,2]*z[:,3]*L1*L2*np.cos(z[:,0]-z[:,1])) -m2*g*(L1*np.cos(z[:,0])+L2*np.cos(z[:,1]))
E = E1 + E2

plt.plot(t, E1, c='b', label='Energía $m_1$')
plt.plot(t, E2, c='r', label='Energía $m_2$')
plt.plot(t, E, c='green', label='Energía sistema')
plt.legend(loc='best')
plt.xlabel('tiempo (s)', fontsize=12)
plt.ylabel('Energía (J)', fontsize=12)
plt.title('Energías en función del tiempo')
plt.show()

#Plano z=0:
X1 = L1*np.sin(z[:,0])
Y1 = -L1*np.cos(z[:,0])
X2 = X1 + L2*np.sin(z[:,1])
Y2 = Y1 - L2*np.cos(z[:,1])

plt.plot(X1, Y1, c='orange', label='Posición de $m_1$', linewidth=0.5)
plt.plot(X2, Y2, c='green', label='Posición de $m_2$', linewidth=0.5)
plt.legend(loc='best')
plt.title('Movimiento de $m_1$ y  $m_2$ en el plano x-y')
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.show()
"""
#Animacion:
Llong=(L1+L2)*1.1

fig, ax3 = plt.subplots()
ax3 = plt.axes(xlim=(-Llong,Llong), ylim=(-Llong,Llong))
ax3.set_xlabel('x (m)')
ax3.set_ylabel('y (m)')

line1,=ax3.plot([],[],lw=2)
line2,=ax3.plot([],[],lw=2)
line3,=ax3.plot([],[],lw=1)
bob1 = plt.Circle((1, 1),Llong*0.02, fc='b')
bob2 = plt.Circle((1, 1),Llong*0.02, fc='r')
time_template = 'time = %.1fs'
time_text = ax3.text(0.05, 0.9, '', transform=ax3.transAxes)

def init():
    bob1.center = (1, 1)
    ax3.add_artist(bob1)
    bob2.center = (0,0)
    ax3.add_artist(bob2)
    line1.set_data([],[])
    line2.set_data([],[]) 
    line3.set_data([],[])
    time_text.set_text('')
    return bob1,bob2,line1,line2,line3,time_text


def animate(i):
    x1, y1 = bob1.center
    x1 = L1*np.sin(z[i,0])
    y1 = -L1*np.cos(z[i,0])
    line1.set_data((0,x1),(0,y1))
    bob1.center = (x1, y1)
    x2, y2 = bob2.center
    x2 = x1+L2*np.sin(z[i,1])
    y2 = y1-L2*np.cos(z[i,1])
    line2.set_data((x1,x2),(y1,y2))
    line3.set_data(L1*np.sin(z[0:i,0])+L2*np.sin(z[0:i,1]),-L1*np.cos(z[0:i,0])-L2*np.cos(z[0:i,1]))
    bob2.center = (x2, y2)
    time_text.set_text(time_template%(i*dt))
    return bob1,bob2,time_text

anim = animation.FuncAnimation(fig, animate, 
                               init_func=init, 
                               frames=10000,
                               interval=5)


plt.show()

"""
########################################################################
#         ANGULOS EN FUNCIÓN DEL TIEMPO PARA EL CAOS                   #
########################################################################

# Ángulos iniciales en grados
theta1_0_deg=20
theta2_0_deg=20
theta3_0_deg=20.00001
theta4_0_deg=20.00001

# Angulos iniciales en radianes
theta1_0=theta1_0_deg*np.pi/180.0
theta2_0=theta2_0_deg*np.pi/180.0
theta3_0=theta3_0_deg*np.pi/180.0
theta4_0=theta4_0_deg*np.pi/180.0

# z0 = [θ10, θ20, ω10, ω20] # Consideremos que no tienen ω0
z0 = [theta1_0,theta2_0,0.0,0.0]
z1 = [theta3_0,theta4_0,0.0,0.0]
t=np.linspace(0,tf,nt)
abserr = 1.0e-8
relerr = 1.0e-6
z=odeint(double_pendulum,z0,t,args=(par,),atol=abserr, rtol=relerr)
z2=odeint(double_pendulum,z1,t,args=(par,),atol=abserr, rtol=relerr)
plt.close('all')


ax1.set_xlabel('tiempo (s)')
ax1.set_ylabel('$\\theta_1$(rad)', color='orange', fontsize=15)
ax2.set_ylabel('$\\theta_2$(rad)', color='blue', fontsize=15)
plt.plot(t, t*0, '--', c='black')
ax1.set_xlim(xmin=0.,xmax=10) #limites del eje x
plt.plot(t[:],z[:,0], linewidth=2, color='blue', label='$\\theta_{1}(0) = 20°$')
plt.plot(t[:],z2[:,0], linewidth=2, color='orange', label='$\\theta_{1}(0) = 20.00001°$')
plt.legend(loc='upper left')
plt.xlabel('Tiempo (s)', fontsize=10)
plt.ylabel('$\\theta_{1}$ (rad)', fontsize=13)

plt.show()


########################################################################
#     CAOS Y ALEJAMIENTO DE LOS PUNTOS EN ÁNGULOS GRANDES              #
########################################################################

theta_inicial = np.linspace(0.0, np.pi, 200)
theta_final = np.linspace(0.0, np.pi, 200)

for i in range(200):
    z0 = [theta_inicial[i], theta_inicial[i], 0.0, 0.0]
    t = np.linspace(0, tf, nt)
    z = odeint(double_pendulum, z0, t, args=(par,), atol=abserr, rtol=relerr)
    theta_final[i] = z[-1, 0]

  
plt.plot(theta_inicial, theta_final, 'o', c='blue', linewidth=0.25)
plt.xlabel('$\\theta_0$ (rad)', fontsize=15)
plt.ylabel('$\\theta_f $ (rad)', fontsize=15)
plt.title('Caos')
plt.show()
"""