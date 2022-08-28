from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as p
L=5
err=0
fig =p.figure(figsize=p.figaspect(1))
N=20
err=1.0
def V(x,y):
	return 0
def ground():
	psi=np.zeros((N,N))
	h=1.0/(N/2.0)
	k=h
	for i in range(1,int(N/2+1)):
		psi[i: -i,i:-i]=k
		k=k+h
	return psi
def st1():
	Z1=np.zeros((int(N/2),N))
	h=1.0/(N/2.0)
	k=h
	for i in range(1,int(N/2+1)):
		Z1[i: -i,i:-i]=k
		k=k+h
	Z2=np.zeros((int(N/2),N))
	h=1.0/(N/2.0)
	k=h
	for i in range(1,int(N/2+1)):
		Z2[i: -i,i:-i]=-k
		k=k+h
	psi=np.concatenate((Z1, Z2), axis=0)
	return psi
def st2():
	Z1=np.zeros((N,int(N/2)))
	h=1.0/(N/2.0)
	k=h
	for i in range(1,int(N/2+1)):
		Z1[i: -i,i:-i]=k
		k=k+h
	Z2=np.zeros((N,int(N/2)))
	h=1.0/(N/2.0)
	k=h
	for i in range(1,int(N/2+1)):
		Z2[i: -i,i:-i]=-k
		k=k+h
	#Z=np.zeroes((N,N))
	psi=np.concatenate((Z1, Z2), axis=1)
	return psi
def sec():
	Z1=np.zeros((int(N/2),int(N/2)))
	h=1.0/(N/2.0)
	k=h
	for i in range(1,int(N/2+1)):
		Z1[i: -i,i:-i]=k
		k=k+h
	Z2=np.zeros((int(N/2),int(N/2)))
	h=1.0/(N/2.0)
	k=h
	for i in range(1,int(N/2+1)):
		Z2[i: -i,i:-i]=-k
		k=k+h
	Z3=np.zeros((int(N/2),int(N/2)))
	h=1.0/(N/2.0)
	k=h
	for i in range(1,int(N/2+1)):
		Z3[i: -i,i:-i]=-k
		k=k+h
	Z4=np.zeros((int(N/2),int(N/2)))
	h=1.0/(N/2.0)
	k=h
	for i in range(1,int(N/2+1)):
		Z4[i: -i,i:-i]=k
		k=k+h
	psi1=np.concatenate((Z1, Z2), axis=1)
	psi2=np.concatenate((Z3, Z4), axis=1)
	psi=np.concatenate((psi1, psi2), axis=0)
	return psi
psi=st1()
E=0
x=y=0
h=L/N
Hpsi=np.zeros((N,N))
x=x+h
y=y+h
for i in range(1,N-1):
	for j in range(1,N-1):
		Hpsi[i,j]=-(psi[i+1,j]+psi[i-1,j]+psi[i,j+1]+psi[i,j-1]-4*psi[i,j])/h**2+V(x,y)*psi[i,j]
		y=y+h
	x=x+h
ExE=0
sum=0
for i in range(1,N-1):
	for j in range(1,N-1):
		sum+=psi[i,j]*psi[i,j]
		ExE+=psi[i,j]*Hpsi[i,j]
ExE=ExE/sum
while(abs(E-ExE)>0.0001 and err>0.0001):
	psic=psi.copy()
	x=0+h
	y=0+h
	for i in range(1,N-1):
		for j in range(1,N-1):
			psi[i,j]=(psi[i+1,j]+psi[i-1,j]+psi[i,j+1]+psi[i,j-1])/(4-h**2*(ExE-V(x,y)))
			x=x+h
		y=y+h
	E=ExE
	Hpsi=np.zeros((N,N))
	x=0+h
	y=0+h
	for i in range(1,N-1):
		for j in range(1,N-1):
			Hpsi[i,j]=-(psi[i+1,j]+psi[i-1,j]+psi[i,j+1]+psi[i,j-1]-4*psi[i,j])/h**2+V(x,y)*psi[i,j]
			x=x+h
		y=y+h
	ExE=0
	sum=0
	for i in range(1,N-1):
		for j in range(1,N-1):
			sum+=psi[i,j]*psi[i,j]
			ExE+=psi[i,j]*Hpsi[i,j]
	ExE=ExE/sum
	err=np.sum(np.abs(psic-psi))
print(E)
psi=psi/(np.sqrt(sum))
x = np.linspace(0,L, N)
y = np.linspace(0,L, N)
X, Y = np.meshgrid(x, y)
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.plot_surface(X, Y, psi,cmap='plasma', edgecolor='none')
p.xlabel("X")
p.ylabel("Y")
p.title("1st Excited Degenrative state1")
ax.set_zlabel("Psi")
psi=st2()
E=0
x=y=0
h=L/N
Hpsi=np.zeros((N,N))
x=x+h
y=y+h
for i in range(1,N-1):
	for j in range(1,N-1):
		Hpsi[i,j]=-(psi[i+1,j]+psi[i-1,j]+psi[i,j+1]+psi[i,j-1]-4*psi[i,j])/h**2+V(x,y)*psi[i,j]
		y=y+h
	x=x+h
ExE=0
sum=0
for i in range(1,N-1):
	for j in range(1,N-1):
		sum+=psi[i,j]*psi[i,j]
		ExE+=psi[i,j]*Hpsi[i,j]
ExE=ExE/sum
while(abs(E-ExE)>0.0001 and err>0.0001):
	psic=psi.copy()
	x=0+h
	y=0+h
	for i in range(1,N-1):
		for j in range(1,N-1):
			psi[i,j]=(psi[i+1,j]+psi[i-1,j]+psi[i,j+1]+psi[i,j-1])/(4-h**2*(ExE-V(x,y)))
			x=x+h
		y=y+h
	E=ExE
	Hpsi=np.zeros((N,N))
	x=0+h
	y=0+h
	for i in range(1,N-1):
		for j in range(1,N-1):
			Hpsi[i,j]=-(psi[i+1,j]+psi[i-1,j]+psi[i,j+1]+psi[i,j-1]-4*psi[i,j])/h**2+V(x,y)*psi[i,j]
			x=x+h
		y=y+h
	ExE=0
	sum=0
	for i in range(1,N-1):
		for j in range(1,N-1):
			sum+=psi[i,j]*psi[i,j]
			ExE+=psi[i,j]*Hpsi[i,j]
	ExE=ExE/sum
	err=np.sum(np.abs(psic-psi))
print(E)
psi=psi/(np.sqrt(sum))
x = np.linspace(0,L, N)
y = np.linspace(0,L, N)
X, Y = np.meshgrid(x, y)
ax = fig.add_subplot(1, 2, 2, projection='3d')
ax.plot_surface(X, Y, psi,cmap='plasma', edgecolor='none')
p.xlabel("X")
p.ylabel("Y")
p.title("1st Excited Degenrative state2")
ax.set_zlabel("Psi")
p.show()



	
	