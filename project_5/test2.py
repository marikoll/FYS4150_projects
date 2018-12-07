#!/usr/bin/python
#Solve vertically propagating Rossby wave problem.
from pylab import *

#Earth's radius
a=6.31e6
#g (m/s**2)
g=9.8
#Reference latitude
phi0 = 45;
#density scale height
hscale = 8.2e3
#reference latitude
phi0 = 45
#mountain wavenumber
k=2*pi/12e6
#top of domain
hd=200.0e3
#damping timescale (seconds, 86400 seconds = 1 day)
tau=20*86400.0
#reference coriolis parameter
f0=2*2*pi*sin(phi0*pi/180)/(86400)
#reference beta
beta=2*2*pi*cos(phi0*pi/180)/(86400*a)
#psi amplitude at surface
c0 = 1e7
#psi amplitude at toa
cnp1 = 0
#number of layers
nlayers=300				
#layer thickness
delta=hd/nlayers
#interface height
z=arange(0,hd+delta,delta)
n = nlayers-1;

def set_arrays(n,delta,tau,k,f0,beta,hscale,c0,cnp1,z):
  a = zeros((n,n),'complex64');
  f = zeros(n,'complex64');
  for jj in range(n):
    zj = z[jj];
    zm = zj-delta/2;

    zp = zj+delta/2;
    #in this code "1j" symbolizes the complex "i"
    #so u(zj)-1j/tau*k corresponds to U-i/(\tau\kappa) in the assignment
    am = (u(zj)-1j/(tau*k))/(rhos(hscale,zj)*delta**2);
    am = am*rhos(hscale,zm)*f0**2/nbv2(zm);

    ap = (u(zj)-1j/(tau*k))/(rhos(hscale,zj)*delta**2);
    ap = ap*rhos(hscale,zp)*f0**2/nbv2(zp);

    aj = -am-ap+ dqdy(delta,f0,beta,hscale,zj)-k**2*u(zj);
    #print jj,zj,zm,zp,am,ap,aj,a[jj,jj]
    a[jj,jj] = aj;
    #print jj,zj,zm,zp,am,ap,aj,a[jj,jj]
    
    if (jj>0):
      a[jj,jj-1]=am;
    else:
      f[jj] = -am*c0
    if (jj<n-1):
      a[jj,jj+1]=ap;
    else:
      f[jj] = -ap*cnp1
  return a,f

#basic state zonal wind
def u(z):
  #u_out =100;
  #u_out = 5+1e-3*z;
  #u_out = max(5,-2+30/cosh((z/1000-10)/10)**2);
  #u_out = max(5,-2+30/cosh((z/1000-30)/20)**2);
  #u_out = max(-1,-2+30/cosh((z/1000-30)/20)**2);
  u_out = 5;
  #u_out=-5;
  #u_out = 90;
  return u_out


#reference density
def rhos(hscale,z):
  rho_out = exp(-z/hscale);
  return rho_out

#square brunt-vaisala frequency
def nbv2(z):
  nbv2_out = 2e-4;
  #n2trop = 2e-4;
  #n2strat = 4e-4;
  #if (z/1000>15):
  #nbv2_out = n2strat
  #else:
  #nbv2_out = 2e-4
  return nbv2_out

#basic state PV
def dqdy(delta,f0,beta,hscale,z):
  dqdy_out = rhos(hscale,z+delta/2)*f0**2/nbv2(z+delta/2)*(u(z+delta)-u(z));
  dqdy_out = dqdy_out - rhos(hscale,z-delta/2)*f0**2/nbv2(z-delta/2)*(u(z)-u(z-delta));
  dqdy_out = dqdy_out/(rhos(hscale,z)*delta**2);
  dqdy_out = beta-dqdy_out;
  return dqdy_out

#m**2 with damping
def m2(z,f0,delta,beta,hscale,tau,k):
  m2_out = nbv2(z)/f0**2*((dqdy(delta,f0,beta,hscale,z)-k**2*u(z))/(u(z)-1j/(tau*k)))-1/4/hscale**2;
  return m2_out

#m**2 without damping
def m2_nodamp(z,f0,delta,beta,hscale,k):
  m2_nodamp_out = nbv2(z)/f0**2*(dqdy(delta,f0,beta,hscale,z)/u(z)-k**2)-1/4/hscale**2;
  return m2_nodamp_out



l,rhs = set_arrays(n,delta,tau,k,f0,beta,hscale,c0,cnp1,z[1:-1])
#solve problem
c=dot(inv(l),rhs)
#m2 = nbv2(1)/f0**2*(beta/u(1)-k**2-(f0/(2*hscale*sqrt(nbv2(1))))**2)
m2_plot = zeros(len(z));
m2_nodamp_plot = zeros(len(z))
u_plot = zeros(len(z))
cgt = zeros(len(z))
for jj in range(len(z)):
  m2_plot[jj] = m2(z[jj],f0,delta,beta,hscale,tau,k);
  u_plot[jj] = u(z[jj]);
  m2_nodamp_plot[jj]= m2_nodamp(z[jj],f0,delta,beta,hscale,k);
  cgt[jj] = 2*u(z[jj])**2*sqrt(m2_nodamp_plot[jj])*k*f0**2*tau/(dqdy(delta,f0,beta,hscale,z[jj])*nbv2(z[jj]));


figure(1)
# a vector of length z that includes the boundary terms
cplot=zeros(2+len(c), 'complex64')
cplot[0]=c0
cplot[1:-1]=c
cplot[-1]= cnp1
clf()

plot(imag(cplot),z/1000,'k')
#You can comment out the following lines after you do the test case: WKB
#solution does not work so well
hold('on')
#plot(real(c0*exp((1j*sqrt(m2_plot)+1/(2*hscale))*z)),z/1000,'r.');
#plot(real(c0*exp((1j*sqrt(m2_nodamp_plot)-1/cgt+1/(2*hscale))*z)),z/1000,'g+');
hold('off')
ylim([0,40]) #the computational domain extends quite
              #high, so it's a good idea to limit it for
               #plotting purposes.
#xlim([-1e8,1e8])
xlabel(r'$\psi$')
ylabel('z')
#subplot(3,1,2)
#plot(u_plot,z)
#ylim([0,40e3])
#subplot(3,1,3)
#plot(m2_nodamp_plot,z)
#ylim([0,40e3])
figure(2)
xmax = 2*pi/k;
x = arange(-xmax,xmax+xmax/100,xmax/100)
psi=zeros((len(x),len(z)))
for jj in range(len(z)):
    psi[:,jj]=real(cplot[jj]*exp(1j*k*x))

cons = arange(5e6,5e7*1.1,5e6)
clf()
contour(x/1000,z/1000,transpose(psi),cons,colors='r');
hold('on')
contour(x/1000,z/1000,transpose(psi),-cons,colors='b');
xlim([x[0]/1000,x[-1]/1000])
ylim([0,40])
xlabel('x')
ylabel('z')
grid('on')

show()