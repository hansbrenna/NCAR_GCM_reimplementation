# -*- coding: utf-8 -*-
"""
Created on Wed May 25 08:58:23 2016

@author: hanbre
An attempt to reeimplement the NCAR GCM from Washington and Williamson 1977
"""
from __future__ import  print_function
import sys
import numpy as np
import pandas as pd
import xarray as xr

def eddy_u():
    fa

def d2lambda(psi,i,j,k,dlambda=dlambda):
    return (psi[i+1,j,k]-psi[i-1,j,k])/(2.*dlambda)
    
def d4lambda(psi,i,j,k,dlambda=dlambda):
    return (psi[i+2,j,k]-psi[i-2,j,k])/(4.*dlambda)

def m2lambda(psi,i,j,k):
    return 0.5*(psi[i+1,j,k]+psi[i-1,j,k])

def m4lambda(psi,i,j,k):
    return 0.5*(psi[i+2,j,k]+psi[i-2,j,k])
    
def d2phi(psi,i,j,k,dphi=dphi):
    return (psi[i,j+1,k]-psi[i,j-1,k])/(2.*dphi)
    
def d4phi(psi,i,j,k,dphi=dphi):
    return (psi[i,j+2,k]-psi[i,j-2,k])/(4.*dphi)
    
def m2phi(psi,i,j,k):
    return 0.5*(psi[i,j+1,k]+psi[i,j-1,k])

def m4phi(psi,i,j,k):
    return 0.5*(psi[i,j+2,k]+psi[i,j-2,k])

def d2s(psi,i,j,k,ds=ds):
    return (psi[i,j,k+1]-psi[i,j,k-1])/float(ds)

def m2s(psi,i,j,k):
    return 0.5*(psi[i,j,k+1]+psi[i,j,k-1])

def update_velocities(u,up,uf,v,vp,vf,a,T,R,p,w,z,dt,dlambda,dphi,ds,omega,lon,rho,Fl,Fph,g):
    #update u: (3.57)   
    for k in u.s_stag    :
        for i in u.lambd:
            for j in u.phi:
                  term1 = up[i,j,k] 
                  term2 = (u[i,j,k]/a*np.cos(lon[j]))*(4/3*d2lambda(u,i,j,k,dlambda)-1/3*d4lambda(u,i,j,k,dlambda)) 
                  term3 = (v[i,j,k]/a)*(4/3*d2phi(u,i,j,k,dphi)-1/3*d4phi(u,i,j,k,dphi))
                  term4_temp = w[i,j,k]*m2s(u,i,j,k)
                  term4 = d2s(term4_temp,i,j,k,ds)
                  term5 = u[i,j,k]*d2s(w,i,j,k,ds)
                  term6 = 2*omega*np.sin(lat[j])+u[i,j,k]/a*np.tan(lat[j])*v[i,j,k]
                  term7 = R/(a*np.cos(phi[j]))*(4/3*m2lambda(T,i,j,k)*d2lambda(np.log(p),i,j,k,dlambda)-
                  1/3*m4lambda(T,i,j,k)*d4lambda(np.log(p),i,j,k,dlambda))
                  term8 = g/(a*np.cos(lon[j]))*(4/3*d2lambda(z,i,j,k,dlambda)-1/3*d4lambda(z,i,j,k,dlambda))
                  term9 = 1/rho*Fl
                  uf[i,j,k]=term1+2*dt*(-term2-term3-term4+term5+term6-term7-term8+term9)
                  
    #update v: (3.58)
    for k in v.s_stag:
        for i in v.lambd:
            for j in v.phi:
                  term1 = vp[i,j,k] 
                  term2 = (u[i,j,k]/a*np.cos(lon[j]))*(4/3*d2lambda(v,i,j,k,dlambda)-1/3*d4lambda(v,i,j,k,dlambda)) 
                  term3 = (v[i,j,k]/a)*(4/3*d2phi(v,i,j,k,dphi)-1/3*d4phi(v,i,j,k,dphi))
                  term4_temp = w[i,j,k]*m2s(v,i,j,k)
                  term4 = d2s(term4_temp,i,j,k,ds)
                  term5 = 4[i,j,k]*d2s(w,i,j,k,ds)
                  term6 = 2*omega*np.sin(lat[j])+u[i,j,k]/a*np.tan(lat[j])*u[i,j,k]
                  term7 = R/a*(4/3*m2phi(T,i,j,k)*d2phi(np.log(p),i,j,k,dphi)-
                  1/3*m4phi(T,i,j,k)*d4phi(np.log(p),i,j,k,dphi))
                  term8 = g/a*(4/3*d2phi(z,i,j,k,dphi)-1/3*d4phi(z,i,j,k,dphi))
                  term9 = 1/rho*Fph
                  vf[i,j,k]=term1+2*dt*(-term2-term3-term4+term5+term6-term7-term8+term9)
    return uf,u,up
    
def pressure_tendency(pf,p,pp,dpT,g,w,a,lon,u,v,rhos,ds,dlambda,dphi,dt):
    rhosu = rhos*u  
    rhosv = rhos*v
    for k in p.s:
        for i in p.lambd:
            for j in p.phi:
                term1 = dpT[i,j]
                term2 = g*m2s(rhos,i,j,k)*w[i,j,k]
                term3 = 0.0                
                for l in xrange(k+1,u.s_stag):
                    term3 += 1/(a*np.cos(lon[j]))*((4/3*d2lambda(rhosu,i,j,l) -
                    1/3*d4lambda(rhosu,i,j,l)) + 
                    (4/3*d2phi(rhosv,i,j,l)-1/3*d4phi(rhosv,i,j,l))*np.cos(lon[j]))
                term3 *= (g*ds)
                pf[i,j,k] = pp[i,j,k]+2*dt*(term1+term2-term3)
                
def prog_humidity(qf,q,qp,u,v,w,a,rho,M,Ep,lon,ds,dlambda,dphi,dt):
    m2sq = q.copy(deep=True)    
    for k in q.s_stag:
        for i in q.lambd:
            for j in q.phi:
                m2sq[i,j,k] = m2s(q,i,j,k)
                
    for k in q.s_stag:
        for i in q.lambd:
            for j in q.phi:
                term1 = u[i,j,k]/(a*np.cos[j])*(4/3.*d2lambda(q,i,j,k)-1/3.*d4lambda(q,i,j,k))
                term2 = v[i,j,k]/a*(4/3.*d2phi(q,i,j,k)-1/3.*d4phi(q,i,j,k))
                term3 = d2s(w*m2sq,i,j,k)
                term4 = q[i,j,k]*d2s(w,i,j,k)
                qf = qp + 2*dt*(-term1-term2-term3+term4+1/rho*M+1/rho*Ep)
    return qf
    
def diag_J(J,u,v,q,p,g,rhos,a,ds,dlambda,dphi):
    for k in J.s:
        term1 = (u.transpose()/(a*np.cos(lon)).transpose()
        term2 = 4/3.*
        
def d2l(psi):
    psi.values=psi.roll

#Constants
pi = np.pi
a = 1
#initializations

nlambda = 73 #even
nphi = 34 #even
dlambda = 2*pi/nlambda
dphi = pi/nphi
lambd = np.arange(0,nlambda)
phi = np.arange(0,nphi)
lon = np.linspace(-pi,pi,nlambda)
lat = np.linspace(-pi/2+0.5*dphi,pi/2-0.5*dphi,nphi)
nlayers = 5
ds = 1./nlayers
s = np.arange(0,nlayers)
s_stag = np.arange(0,nlayers)
lev = np.linspace(0,1,nlayers)
lev_stag = np.zeros([nlayers])
lev_stag[1:]=np.linspace(0.5*ds,1-0.5*ds,nlayers-1)

u_n = xr.DataArray(data=np.zeros([nlambda,nphi,nlayers]),dims=('lambd','phi','s_stag'),coords={'lambd':lambd,'phi':phi,'s_stag':s_stag})
u_nm1 = u_n.copy(deep=True)
u_np1 = u_n.copy(deep=True)