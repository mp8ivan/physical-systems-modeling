import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import matplotlib.animation as animation


fig=plt.figure()
fig.set_dpi(100)

#Datos a introducir:
L1=0.1      #longitud del péndulo 1
L2=0.1      #longitud del pendulo 2
m1=1.0      #masa del pendulo 1
m2=1.0      #masa del pendulo 2
g=9.81      #aceleracion de la gravedad
tf=10.0     #tiempo de simulacion

m12=m1+m2
par=[L1,L2,m1,m2,m12,g,0.2,0.2]

#Ecuaciones del movimiento:
def damped_double_pendulum(z,t,par): #devuelve theta1, theta2, thetapunto1, thetapunto2 en columnas
    z1,z2,z3,z4=z  
    
    L1,L2,m1,m2,m12,g, k1, k2 = par

    sinz=np.sin(z1-z2)
    cosz=np.cos(z1-z2)
    
    sin2z = np.sin(2*(z1-z2))
    
    sinz1=np.sin(z1)
    cosz2 = np.cos(z2)
    
    z42=z4*z4
    z32=z3*z3
    
    gamma1 = 2*k1*z3- 2*k2*z4*cosz
    gamma2 = 2*k1*z3*cosz - 2*(m1+m2)/m2*k2*z4

    dzdt=[z3,z4,(m2*L1*z32*sin2z + 2*m2*L2*z42*sinz + 2*g*m2*cosz2*sinz + 2*g*m1*sinz1 + gamma1)/(-2*L1*(m1+m2*sinz*sinz)),
          (m2*L2*z42*sin2z + 2*(m1+m2)*L1*z32*sinz + 2*g*(m1+m2)*cosz2*sinz + gamma2)/(2*L2*(m1+m2*sinz*sinz))]
    return dzdt


#Llamada a odeint que resuelve las ecuaciones de movimiento:
nt=500  #numero de intervalos de tiempo
dt=tf/nt


#Angulos iniciales en grados:
theta1_0_deg = 20.0
theta2_0_deg = 20.0

#Angulos iniciales en radianes:
theta1_0=theta1_0_deg*np.pi/180.0
theta2_0=theta2_0_deg*np.pi/180.0
z0=[theta1_0,theta2_0,0.0,0.0] # Valores iniciales

t=np.linspace(0,tf,nt)

abserr = 1.0e-8
relerr = 1.0e-6
z=odeint(damped_double_pendulum,z0,t,args=(par,),atol=abserr, rtol=relerr)
plt.close('all')

#Angulos en funcion del tiempo:
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.set_xlabel('Tiempo (s)')
ax1.set_ylabel('$\\theta_1$ (rad)', color='blue', fontsize=15)
ax2.set_ylabel('$\\theta_2$ (rad)', color='orange', fontsize=15)
ax1.tick_params('y', colors='blue')
ax2.tick_params('y', colors='orange')
ax1.set_xlim(xmin=0.,xmax=tf) #limites del eje x
line1, = ax1.plot(t[:],z[:,0], linewidth=2, color='blue')
line2, = ax1.plot(t[:],z[:,1], linewidth=2, color='orange')
plt.show()

#Energia de cada masa y del sistema:
K1 = 0.5*m1*z[:,2]**2*L1**2
U1 = -m1*g*L1*np.cos(z[:,0])
E1 = 0.5*m1*z[:,2]**2*L1**2 - m1*g*L1*np.cos(z[:,0])
K2 = 0.5*m2*(z[:,2]**2*L1**2 + z[:,3]**2*L2**2 + 2*z[:,2]*z[:,3]*L1*L2*np.cos(z[:,0]-z[:,1]))
U2 = -m2*g*(L1*np.cos(z[:,0])+L2*np.cos(z[:,1]))
E2 = 0.5*m2*(z[:,2]**2*L1**2 + z[:,3]**2*L2**2 + 2*z[:,2]*z[:,3]*L1*L2*np.cos(z[:,0]-z[:,1])) - m2*g*(L1*np.cos(z[:,0])+L2*np.cos(z[:,1]))
E = E1 + E2

plt.plot(t, K2, c='blue', label='Energía cinética $m_2$')
plt.plot(t, U2, c='orange', label='Energía potencial $m_2$')
plt.plot(t, E2, c='yellow', label='Energía total $m_2$')
plt.legend(loc='best')
plt.xlabel('Tiempo (s)', fontsize=12)
plt.ylabel('Energía (J)', fontsize=12)
plt.show()

#Animación:
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

#Caos y alejamiento de los puntos en ángulos grandes:
theta_inicial = np.linspace(0.0, np.pi, 200)
theta_final = np.linspace(0.0, np.pi, 200)

for i in range(200):
    z0 = [theta_inicial[i], theta_inicial[i], 0.0, 0.0]
    t = np.linspace(0, tf, nt)
    z = odeint(damped_double_pendulum, z0, t, args=(par,), atol=abserr, rtol=relerr)
    theta_final[i] = z[-1, 0]

  
plt.plot(theta_inicial, theta_final, 'o', c='blue', linewidth=0.01)
plt.axhline(0, linestyle='--', color='black')
plt.xlabel('$\\theta_{10}$ (rad)', fontsize=15)
plt.ylabel('$\\theta_1$ (rad)', fontsize=15)
plt.ylim(-5, 5)
plt.show()