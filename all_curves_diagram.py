#!/usr/bin/env python
# FITTING
from pylab import *
from scipy.optimize import curve_fit
from scipy.integrate import quad

#params = {'legend.fontsize': 12, #}
#          'xtick.minor.size': 6,'ytick.minor.size': 6}
fig_size = [700/72.27 ,520/72.27]
#fig_size = [520/72.27 ,520/72.27]
params = {'axes.labelsize': 22, 'legend.fontsize': 16,
          'xtick.labelsize': 22, 'ytick.labelsize': 22,
          'xtick.major.size': 14,'ytick.major.size': 14,
          'xtick.minor.size': 7,'ytick.minor.size': 7,
          'xtick.direction': 'in', 'ytick.direction': 'in',
          'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
          'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0,
          'text.usetex': False, 'figure.figsize': fig_size, 'axes.linewidth': 2,
          'xtick.major.pad': 5,
          'ytick.major.pad': 10,
          'figure.subplot.bottom': 0.110,'figure.subplot.top': 0.975,'figure.subplot.left': 0.140,'figure.subplot.right': 0.977}

rcParams.update(params)

def sec(x): return 1.0/(cos(x))

def rCircle(phi, R):
 r_hat = array( [sin(phi), 1-cos(phi)] )
 return R*r_hat

def rCycloid(phi, r):
 r_hat = array( [ 2*phi + sin(2*phi), 1-cos(2*phi)] )
 return r*r_hat

def rCatenary(phi, C):
 alpha = arcsinh( tan(phi) )
 r_hat = array( [ alpha, cosh(alpha)-1 ] )
 return C*r_hat

def rParabola(phi, P):
 r_hat = array( [tan(phi)/2, (tan(phi)/2)**2] )
 return P*r_hat




markers = [  'p', 'o', 'H', 'd',      'D', '<', '^', 's', 'h', 'H', 'X', '*', 'P', 'd', '|']
colors = ["red",
          "#56B4E9",  # sky blue
          "#E69F00",  # golden orange
          "#CC79A7",  # magenta
          "#009E73",  # teal green
          'magenta',
          '#4B0082',  # dark violet
          '#0072B2', # royal blue
          'blue']



fig = figure(1)
ax = subplot(111)
xlabel(r'$x(\varphi)/H$')
ylabel(r'$y(\varphi)/H$')

H=1

r=0.5*H
phi = linspace(0, pi/2)
x,y = rCycloid(phi,r)
ax.plot( x, y, color='k', lw=6, label='Cycloid', zorder=0 ) 
ax.plot(x[-1],y[-1], 'o', mec='k',mew=2, mfc='grey',  ms=15, zorder=7)
phi = linspace(pi/2, 1.5*pi/2)
x,y = rCycloid(phi,r)
ax.plot( x, y, '--', color='k', dashes=[10,3,3,3], lw=1,  zorder=0 ) 


phi_cycloid = arcsin(4*pi/(4+pi*pi))
phi = linspace(0, phi_cycloid)
R=(4+pi*pi)/8 * H
x,y = rCircle(phi,R)
ax.plot( x, y,  color='red', lw=4, label='Circle', zorder=-1) 
phi = linspace(phi_cycloid, pi/2)
x,y = rCircle(phi,R)
ax.plot( x, y,'--', color='red', lw=1, zorder=-1) 



alpha_cycloid =  1.14319
phi = linspace(0, arctan(sinh(alpha_cycloid)))
C=pi/(2*alpha_cycloid) * H
x,y = rCatenary(phi,C)
ax.plot( x, y, color='blue', lw=4, label='Catenary', zorder=-2) 
phi = linspace(arctan(sinh(alpha_cycloid)), 1.4*pi/4)
x,y = rCatenary(phi,C)
ax.plot( x, y, '--', color='blue', dashes=[10,2], lw=1, zorder=-2) 


phi_parabola = arctan(4/pi)
phi = linspace(0, phi_parabola)
P= pi*pi * H/4
x,y = rParabola(phi,P)
ax.plot( x, y,  color='darkviolet', dashes=[5,0.5,0.5,0.5], lw=4, label='Parabola', zorder=-3) 
phi = linspace(phi_parabola, 1.5*phi_parabola)
P= pi*pi * H/4
x,y = rParabola(phi,P)
ax.plot( x, y, '--',  color='darkviolet', dashes=[20,4], lw=1,  zorder=-3) 


ax.set_ylim(0,1.2*H)
#ax.set_ylim(0,1.1*pi/2)
ax.set_xlim(0,1.1*pi/2)
ax.set_xticks([0,0.5,pi/4, 1, pi/2], [0,0.5,'$\pi/4$', 1, '$\pi/2$' ])
ax.set_aspect('equal')

ax.legend()

#savefig('all_curves.png')
show()			
