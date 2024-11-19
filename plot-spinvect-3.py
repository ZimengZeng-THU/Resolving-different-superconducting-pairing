#!~/.conda/envs/nlp/bin/python3.7
import numpy as np
import matplotlib.pyplot as plt
import math as mmm
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import string
import matplotlib.mlab as ml
import matplotlib
from scipy.interpolate import Rbf
import pickle
nk1=20
nk2=20
nk3=20
prefix='Pb'

bxsffile=prefix+'.bxsf'
Vkkfile= prefix+'.lambda_kkq'
datafile= prefix+'.lambda_FS'



bxsf=open(bxsffile)
#kmesh0=np.loadtxt('kmesh41.dat')
lines = bxsf.readlines()#获取所有行
sum = 0
Ek = []
nl=0
for line in lines:#第i行
    nl=nl+1
    
    if nl >20:
        
        if line[1]!='B' and line[1]!='E' and line[2]!='E':
            line=float(line)
            #print(line)
            Ek.append(line)
print(len(Ek))
N=nk1+1
NN=nk1
xk=0
kmesh0=np.zeros((N*N*N,3))
while xk<=NN:
    yk=0
    while yk<=NN:
        zk=0
        while zk<=NN:
            kmesh0[zk+yk*N+xk*N*N]=[xk/N,yk/N,zk/N]
            zk=zk+1
        yk=yk+1
    xk=xk+1
#print(kmesh0[zk-1+yk*N-N+xk*N*N-20*N*N])
print(kmesh0)
print(kmesh0[0,2])
nbwan=int(len(Ek)/len(kmesh0))
Lk=len(kmesh0)
print(nbwan)
kmesh=np.array(kmesh0)
nn=0
print(np.concatenate((kmesh, np.array(kmesh0)),axis=0))
#print(X1)
while nn<nbwan-1:
    kmesh=np.concatenate((kmesh, np.array(kmesh0)),axis=0)
    nn=nn+1
print(len(kmesh))






Note_to_pauli_mat='''
wfc_type==1: (a1_up,a1_dn, a2_up,a2_dn,   ...)
wfc_type==2: (a1_up,a2_up,... a1_dn,a2_dn,...)
'''

def pauli_mat(norb,wfc_type=1):
    sigma_x=np.array([[0., 1. ], [1,  0.]],complex)
    sigma_y=np.array([[0.,-1.j], [1.j,0.]],complex)
    sigma_z=np.array([[1., 0. ], [0.,-1.]],complex)
    sigma_i=np.array([[1., 0. ], [0., 1.]],complex)
    II=np.diag(np.ones(norb,complex))
    if wfc_type==1:
        pauli_x=np.kron(II,sigma_x)
        pauli_y=np.kron(II,sigma_y)
        pauli_z=np.kron(II,sigma_z)
        pauli_i=np.kron(II,sigma_i)
    elif wfc_type==2:
        pauli_x=np.kron(sigma_x,II)
        pauli_y=np.kron(sigma_y,II)
        pauli_z=np.kron(sigma_z,II)
        pauli_i=np.kron(sigma_i,II)
    #print ( 'wfc type {0}'.format(wfc_type),)
    return pauli_x,pauli_y,pauli_z,pauli_i


def write_spin_texture(calc_bands,spinors,k_vec=None,evals=None):
    nkpt=spinors.shape[1]
    fmt = '{:>10s} '*2+'{:>10s} '*3+'\n'
    with open('spin_texture_tb.dat','w') as fw:
        if k_vec is not None and evals is not None: 
            fw.write(fmt.format('k_dist','band_en','Sx','Sy','Sz'))
            for ib,iband in enumerate(calc_bands):
                for ik in range(nkpt):
                    fw.write(('{:10.3f}'*3+' {:10.3f} ').format(k_vec[ik,0],k_vec[ik,1],k_vec[ik,2],evals[iband,ik]))
                    fw.write(('{:10.5f} '*3+'\n').format(*tuple(spinors[ib,ik])))
                fw.write('\n')

        else:
            fw.write('# kpt/band index starts from 0\n')
            fw.write(fmt.format('kpt','band','Sx','Sy','Sz'))
            for ib,iband in enumerate(calc_bands):
                for ik in range(nkpt):
                    fw.write('{:10d} {:10d} '.format(ik,iband))
                    fw.write(('{:10.5f} '*3+'\n').format(*tuple(spinors[ib,ik])))
                fw.write('\n')


def calc_spin_texture(evec,wfc_type,calc_bands,write_spin=True,k_vec=None,evals=None):
    shape = evec.shape
    if len(shape)>3: evec=evec.reshape(shape[0],shape[1],shape[2]*shape[3])
    #if np.mod(evec.shape[-1],2)==1: exit( 'Error! Dim of wfc is odd, non spinor case?')
    nband,nkpt,norb=evec.shape
    #if norb//2*2!=norb: Error ('calc_spin_texture: norb must be even!')
    norb=norb//2
    spinors=np.zeros((len(calc_bands),nkpt,3),float)
    pauli=pauli_mat(norb,wfc_type=wfc_type)
    for ib,iband in enumerate(calc_bands):
        for ik in range(nkpt):
            wfc=evec[iband,ik]
            #for ii in range(3): spinors[ib,ik,ii] = np.einsum('j,jk,k',wfc.conj(),pauli[ii],wfc).real
            for ii in range(3): spinors[ib,ik,ii]=np.dot(np.dot(wfc.conj(),pauli[ii]),wfc).real
            norm=np.linalg.norm(spinors[ib,ik])
            if norm<1e-3: print ('ib {0} ik {1} spinors {2}'.format(ib+1,ik+1,spinors[ib,ik]))
            spinors[ib,ik]/=norm*2

    if write_spin: write_spin_texture(calc_bands,spinors,k_vec=k_vec,evals=evals)
    return spinors

def calc_spin_wfc(wfc):
    norb=len(wfc)
    #if norb//2*2!=norb: Error ('calc_spin_texture: norb must be even!')
    norb=norb//2
    spin=np.zeros((3),float)
    pauli=pauli_mat(norb,wfc_type=wfc_type)
    for ii in range(3): spin[ii]=np.dot(np.dot(wfc.conj(),pauli[ii]),wfc).real
    norm=np.linalg.norm(spin)
    #if norm<1e-3: print ('spinors {2}'.format(spin))
    spin/=norm*2
    return spin

def plot_bands(k_dist,k_node,k_label,evals,band_color='C0',band_alpha=1,elim=None,save=False,show=True):
    fig, ax = plt.subplots(figsize=(5,4))
    for i in range(evals.shape[0]):
        ax.plot(k_dist,evals[i],"-",c=band_color,alpha=band_alpha)
    for n in range(len(k_node)):
        ax.axvline(x=k_node[n],linewidth=0.5, color='k')
    ax.set_ylabel("$E$ (eV)")
    ax.set_xlim(k_dist[0],k_dist[-1])
    ax.set_xticks(k_node)
    #ax.set_xticklabels(['$\mathrm{'+item+'}$' for item in k_label])
    if elim is not None: ax.set_ylim(elim)
    fig.tight_layout()
    if save: fig.savefig("TB_band_quick",dpi=300)
    if show: plt.show()
    return fig










#bxsf=np.loadtxt('wannier90.bxsf')

Vkkdata=np.loadtxt(datafile)
LVkk=len(Vkkdata)
kkkk=np.unique(Vkkdata[:,3])
Lk=len(kkkk)
aaaa=np.unique(Vkkdata[:,8])
Libnd=len(aaaa)
Lkk=np.zeros((Lk))



lambdaFS=np.loadtxt(datafile)
kindex=np.unique(lambdaFS[:,3],return_index=True)
kpoints = lambdaFS[kindex[1],0:3]  # This is all the unique x-points
    



Vkkdata1=np.zeros((Lk,Lk,Libnd,Libnd,3),dtype=np.complex64)
Vkkdata2=np.zeros((Lk,Lk,Libnd,Libnd,Libnd,Libnd),dtype=np.complex64)


dosef0=np.zeros((Lk,Libnd))
dosef01=np.zeros((Lk,Libnd))
dosef1=np.zeros((Lk,Libnd))

Vkk=np.zeros((Lk*Libnd*Libnd,Lk*Libnd*Libnd),dtype=np.complex64)
Vkkss=np.zeros((Lk*Libnd*Libnd,Lk*Libnd*Libnd),dtype=np.complex64)

dkk=np.zeros((Lk,Libnd,Libnd),dtype=np.complex64)
dkk0=np.zeros((Lk,Libnd,Libnd),dtype=np.complex64)
gkk=np.zeros((Lk,Libnd,Libnd),dtype=np.complex64)
orderk1=np.zeros(LVkk)
spinor=np.zeros((nbwan,nk1,nk2,nk3,2,2),dtype=np.complex64)
spinorl=np.zeros((nbwan,nk1,nk2,nk3,3),dtype=np.float64)
spinorll=np.zeros((nbwan,nk1,nk2,nk3,3),dtype=np.float64)
#orderk2=np.zeros(LVkk)
phase=np.zeros((nbwan,nk1,nk2,nk3),dtype=np.complex64)
spinorl=np.zeros((nbwan,nk1,nk2,nk3,3),dtype=np.float64)
Ef=11.2876
ne=4
s1=np.zeros((2,2),dtype=np.complex64)
s2=np.zeros((2,2),dtype=np.complex64)
s3=np.zeros((2,2),dtype=np.complex64)
s1[0,0]=0
s1[0,1]=1
s1[1,0]=1
s1[1,1]=0
s2[0,0]=0
s2[0,1]=-1j
s2[1,0]=1j
s2[1,1]=0
s3[0,0]=1
s3[0,1]=0
s3[1,0]=0
s3[1,1]=-1  

wfc_type=1
d1=np.loadtxt('wfc333.dat')
Lk00=nk1*nk2*nk3
#nbwan=int(len(d1)/Lk)
#lwfc=int((len(d1[0,:])-4)/2)
wfck=np.zeros((nbwan,Lk,nbwan),dtype=np.complex64)
#print(lwfc)
LLwfc=len(d1)

for i in range(LLwfc):
    ik=int(d1[i,0]-1)
    ibnd=int(d1[i,1]-1)
    jbnd=int(d1[i,2]-1)
    imode=int(d1[i,3]-1)
    if imode==0:
        wfck[ibnd,ik,jbnd]=d1[i,4]+1j*d1[i,5]

wfc333=wfck[0,0]
wfc3333=wfck[1,0]
www=[wfc333,wfc3333]
print(calc_spin_wfc(wfc333))
diakk=[[1,1],[-1,1]]
newspin=np.dot(np.conjugate(np.transpose(diakk)),www)
newwfc0=newspin[0,:]
newwfc1=newspin[1,:]
print(calc_spin_wfc(newwfc0))
print(calc_spin_wfc(newwfc1))
print('nbwan')
print(nbwan)


def plotdddd(channel):

    eevvv = np.loadtxt('delta'+'/'+'dkk'+str(ne)+'.txt').view(complex).reshape(-1)
    lll=abs(np.sqrt((eevvv**2).sum(axis=0)))
    for i in range(Lk):
        for iii in range(int(Libnd/2)):
            for iibnd in range(2):
                for jjbnd in range(2):
                    dkk0[i][2*iii+iibnd, 2*iii+jjbnd] = eevvv[i * Libnd * Libnd +iii*4+ iibnd * 2 + jjbnd]
        
        #print(dkk[i])


    eevvv_real, eevvv_imag = np.loadtxt('delta'+'/'+'dkk'+str(ne)+'.txt', unpack=True)
    eevvv = eevvv_real + 1j * eevvv_imag
    norm=matplotlib.colors.Normalize(vmin=0,vmax=1)
    eevvv = np.loadtxt('delta'+'/'+'gkk'+str(ne)+'.txt').view(complex).reshape(-1)
    for i in range(Lk):
        for iii in range(int(Libnd/2)):
            for iibnd in range(2):
                for jjbnd in range(2):
                   gkk[i][2*iii+iibnd,2*iii+jjbnd]=eevvv[i*Libnd*Libnd+iii*4+2*iibnd+jjbnd]
        #print(i)
        #print(gkk[i])
        

        

    aa0=0
    aa00=0
    aa1=0
    kline=np.linspace(0,len(dkk[:,0,0]),len(dkk[:,0,0]))
    #print(len(kline))
    #print(len(dkk[:,0,0]))

    finek3=0
    ik=finek3
    
    dkkxyz=np.zeros((nbwan,nk1,nk2,nk3))
    Ekkxyz=np.zeros((nbwan,nk1,nk2,nk3))
    #print(Ekkxyz)
    for ibnd in range(int(nbwan)):
        for xk in range(NN):
            for yk in range(NN):
                for zk in range(NN):
                    index=ibnd*N*N*N+xk*N*N+yk*N+zk
                    #print(Ek[index])
                    Ekkxyz[ibnd,xk,yk,zk]=Ek[index]-Ef
    finek3=0
    ik=finek3
    #lll=0
    cmax=max(abs(dkk0[:,0,0]))
    for i in range(Lk):
        #print(i)
        ik=kindex[1][i]
        kk1=int(lambdaFS[ik,3])-1
        kkk1=int(lambdaFS[ik,4])-1
        #print(kk1)
        #print(dkk[kk1])
        #print(dkk[kkk1])
        #print(abs(gkk[i]))
        
        
        kx=lambdaFS[kindex[1][i],0]
        ky=lambdaFS[kindex[1][i],1]
        kz=lambdaFS[kindex[1][i],2]
        
        nkx=int(kx*nk1)
        nky=int(ky*nk2)
        nkz=int(kz*nk3)
        
        xk=nkx
        yk=nky
        zk=nkz
        cmax=max(abs(dkk0[:,0,0]))
        finek3=nkx*nk2*nk3+nky*nk3+nkz
        sum=0
        for ii in range(int(Libnd/2)):# print([[dkk[kk1][2*ibnd,2*ibnd],dkk[kk1][2*ibnd+1,2*ibnd]],[dkk[kk1][2*ibnd,2*ibnd+1],dkk[kk1][2*ibnd+1,2*ibnd+1]]])
            #print([[dkk[kkk1][2*ibnd,2*ibnd],dkk[kkk1][2*ibnd+1,2*ibnd]],[dkk[kkk1][2*ibnd,2*ibnd+1],dkk[kkk1][2*ibnd+1,2*ibnd+1]]])
            gkkt=np.array([[gkk[kkk1][2*ii,2*ii],gkk[kkk1][2*ii+1,2*ii]],[gkk[kkk1][2*ii,2*ii+1],gkk[kkk1][2*ii+1,2*ii+1]]])
            #print(abs(np.array(gkkt)))
            #if abs(gkkt[0,1])<abs(gkkt[0,0]):
            #     dkk[i]=dkk[i]*-1
            #if gkk[i][]
            #delta=np.real(dkk[i,ii,ii])/cmax
            EEE=lambdaFS[kindex[1][i],9]
   
            
            summ=0
            for iii in range(nbwan):
                index=ii*N*N*N+xk*N*N+yk*N+zk
                #print(kx*nk1,nkx)
                #print(EEE,Ekkxyz[ii,nkx,nky,nkz])
                #有些互为空间反演的k点，只有一个点有非零的gap，这是因为拟合的能带误差，导致把算出来的gap和能带匹配起来时没匹配到
                if abs(Ekkxyz[iii,nkx,nky,nkz]-EEE)<0.002 and iii%2==0:
                    summ=summ+abs(dkk0[kk1][2*ii, 2*ii])**2+abs(dkk0[kk1][2*ii, 2*ii+1])**2+abs(dkk0[kk1][2*ii+1, 2*ii])**2+abs(dkk0[kk1][2*ii+1, 2*ii+1])**2
                    nkxx=nkx
                    nkyy=nky
                    nkzz=nkz
                    if nkz!=0:
                        nkzz=nk1-nkz
                    if nky!=0:
                        nkyy=nk2-nky
                    if nkx!=0:
                        nkxx=nk1-nkx
                    #print(nkxx,nkyy,nkzz)
                    spin=spinor[iii,nkx,nky,nkz]
                    spinkkk=spinor[iii,nkxx,nkyy,nkzz]
                    
                    phase0=phase[iii+1,nkx,nky,nkz]
                    phase0kkk=phase[iii+1,nkxx,nkyy,nkzz]
                    #phase0=np.conjugate(phase0)
                    #phase0kkk=np.conjugate(phase0kkk)
                    phase0=1
                    phase0kkk=1
                    
                    
                    
                    
                    
                    
                    deltakk = np.array([[dkk0[kk1][2*ii, 2*ii], dkk0[kk1][2*ii, 2*ii+1]], [dkk0[kk1][2*ii+1, 2*ii], dkk0[kk1][2*ii+1, 2*ii+1]]])
                    deltakk =deltakk+np.conjugate(np.transpose(deltakk))
                    #deltakk = np.array([[dkk0[kk1][2*ii, 2*ii], dkk0[kk1][2*ii, 2*ii+1]], [np.conjugate(dkk0[kk1][2*ii, 2*ii+1]), np.conjugate(dkk0[kk1][2*ii, 2*ii])]])

                    [eigkk, diakk] = np.linalg.eig(deltakk)
                    #diakk=1
                    newdeltakk =  np.dot(np.dot(np.conjugate(np.transpose(diakk)), deltakk), diakk)
                    
                    wfckk00=wfck[iii,kk1]
                    wfckk01=wfck[iii+1,kk1]
                    #if np.real(wfckk00[0])<0:
                    #    wfckk00=-wfckk00
                    #if np.real(wfckk01[0])<0:
                    #    wfckk01=-wfckk01
                    
                    deltakkk = np.array([[dkk0[kkk1][2*ii, 2*ii], dkk0[kkk1][2*ii, 2*ii+1]], [dkk0[kkk1][2*ii+1, 2*ii], dkk0[kkk1][2*ii+1, 2*ii+1]]])
                    deltakkk =deltakkk+np.conjugate(np.transpose(deltakkk))
                    #deltakkk = np.array([[dkk0[kkk1][2*ii, 2*ii], dkk0[kkk1][2*ii, 2*ii+1]], [np.conjugate(dkk0[kkk1][2*ii, 2*ii+1]), dkk0[kkk1][2*ii+1, 2*ii+1]]])

                    #deltakkk = np.array([[deltakkk[0,0], 1/phase0kkk*deltakkk[0,1]], [1/np.conjugate(phase0kkk)*deltakkk[1,0], deltakkk[1,1]]])

                    [eigkkk, diakkk] = np.linalg.eig(deltakkk)
                    #diakk=1
                    newdeltakkk =  np.dot(np.dot(np.conjugate(np.transpose(diakkk)), deltakkk), diakkk)
                    wfckkk00=wfck[iii,kkk1]
                    wfckkk01=wfck[iii+1,kkk1]
                    #if np.real(wfckk00[0])<0:
                    #    wfckk00=-wfckk00
                    #if np.real(wfckkk01[0])<0:
                    #    wfckkk01=-wfckkk01

                    spinkk=np.dot(np.dot(np.conjugate(np.transpose(diakk)), spinor[iii,nkx,nky,nkz]), diakk)
                    #spin3=spin
                    
                    wfckk=[wfckk00,wfckk01]
                    wfckk=np.dot(np.conjugate(np.transpose(diakk)),wfckk)
                    wfckk0=wfckk[0,:]
                    wfckk1=wfckk[1,:]
                    
                    wfckkk=[wfckkk00,wfckkk01]
                    wfckkk=np.dot(np.conjugate(np.transpose(diakkk)),wfckkk)
                    wfckkk0=wfckkk[0,:]
                    wfckkk1=wfckkk[1,:]

                    skk=calc_spin_wfc(wfckk0)
                    skkk=calc_spin_wfc(wfckkk0)

                    sx=skk[0]
                    sy=skk[1]
                    sz=skk[2]
                    sx1=skkk[0]
                    sy1=skkk[1]
                    sz1=skkk[2]
                    spink=[sx,sy,sz]
                    spinkkk=[sx,sy,sz] 
                    dkk[kk1]=newdeltakk
                    #if sz<0:
                    #    dkk[kk1]=-newdeltakk
                    sum=sum+abs(newdeltakk[0,0])*abs(newdeltakk[0,0])
                    #spin3=spin
                    #[eigsk, diask] = np.linalg.eig(spin3)
                    #newdeltakk =  np.dot(np.dot(np.conjugate(np.transpose(diask)), deltakk), diask)
                    
                    if kx==0.15 and ky==0.0 and kz==0.15:
                        #overlap=abs(np.dot(np.transpose(np.conjugate(wfckk0)),wfckk1))
                        #sx=np.real(spin3[0,1]+spin3[1,0])
                        #sy=np.real(1j*spin3[0,1]-1j*spin3[1,0])
                        #sz=np.real(spin3[0,0]-spin3[1,1])
                        kx=round(kx,3)
                        ky=round(ky,3)
                        print(kx,ky,kz)
                        print(newdeltakk/cmax)
                        print(newdeltakkk/cmax)
                    
                        print(skkk)
                        print('\n')
                    if kx==0.0 and ky==0.15 and kz==0.15:
                        print(kx)
                        print(ky)
                        print(kz)
                        print(newdeltakk/cmax)
                        print(newdeltakkk/cmax)
                        print(skk)
                        print(skkk)
                        print('\n')
                    if iii%2==0:
                        #print(newdeltakk)
                        dkkxyz[iii,nkx,nky,nkz]=np.real(newdeltakk[0,0])
                        if sx<0:
                            dkkxyz[iii,nkx,nky,nkz]=-np.real(newdeltakk[0,0])

                        break
        #print(dkkxyz[ii,xk,yk,zk],dkkxyz[ii,nk1-xk-1,nk2-yk-1,nk3-zk-1])#print(dkkxyz[ii,xk,yk,zk],dkkxyz[ii,nk1-xk-1,nk2-yk-1,nk3-zk-1])
            
            #Ekkxyz[ibnd,nkx,nky,nkz]=EEE
        EEE=lambdaFS[kindex[1][0],9]
    

    
    #print(abs(np.sqrt((eevvv**2).sum(axis=0))))
    print(sum)
    print(summ)
    cmax=max(abs(dkk0[:,0,0]))
    norm=matplotlib.colors.Normalize(vmin=-0.6,vmax=0.6)
    a=0
    for i in range(Lk):
        kx=lambdaFS[kindex[1][i],0]
        ky=lambdaFS[kindex[1][i],1]
        kz=lambdaFS[kindex[1][i],2]
        nkx=int(nk1*kx)
        nky=int(nk2*ky)
        nkz=int(nk3*kz)
        #if kx>0.5:
        #    kx=kx-1
        #if ky>0.5:
        #    ky=ky-1
        if kz>0.5:
            kz=kz-1
        if kz==0:
            b=np.real(dkk[i,0,0])/cmax
            a=np.real(dkk[i,0,0])
            ccc=spinorl[iii,nkx,nky,nkz][1]
            #plt.scatter(kx,ky, c=b,cmap='cool',norm=norm)
            plt.scatter(kx,ky, c=b,cmap='cool',norm=norm)
            #print(kx,ky,np.real(dkk[i,0,0])/cmax)

    plt.colorbar()
    #plt.title('Energy Spectrum of Chain with {} Sites'.format(Nsites))
    plt.show()
    plt.savefig('p.png')
    
    for i in range(nbwan):
        rawdata=[]
        #points=np.array([[0,0,0]])
        x=[]
        y=[]
        z=[]
        
        xf=[]
        yf=[]
        zf=[]
        for nkx in range(nk1):
            for nky in range(nk2):
                for nkz in range(nk3):
                    xf.append(nkx)
                    yf.append(nky)
                    zf.append(nkz)
                    if dkkxyz[i,nkx,nky,nkz]>0.00001 or dkkxyz[i,nkx,nky,nkz]<-0.00001:
                        rawdata.append(dkkxyz[i,nkx,nky,nkz])
                        x.append(nkx)
                        y.append(nky)
                        z.append(nkz)
                        #points=np.concatenate((points,np.array([[nkx,nky,nkz]])),axis=0)
        #print(dkkxyz.shape)
        #print(points.shape)
        points=np.array(np.transpose([x,y,z]))
        mesh=np.array(np.transpose([xf,yf,zf]))
        #print(points)
        if len(points)>=1:
            #grid_z = griddata(points,np.array(rawdata), mesh, method='linear',fill_value='extrapolate')
            fun = Rbf(x,y,z,rawdata,method='linear',smooth=0)
            
            for nkx in range(nk1):
                for nky in range(nk2):
                    for nkz in range(nk3):
                        #print(dkkxyz[i,nkx,nky,nkz],fun(nkx,nky,nkz))
                        dkkxyz[i,nkx,nky,nkz]=fun(nkx,nky,nkz)
                        #dkkxyz[i,nk1-nkx-1,nk2-nky-1,nk3-nkz-1]=-fun(nkx,nky,nkz)
                        ppp=0
                        
    
    for i in range(Lk):

        
        kx=lambdaFS[kindex[1][i],0]
        ky=lambdaFS[kindex[1][i],1]
        kz=lambdaFS[kindex[1][i],2]
        
        nkx=int(kx*nk1)
        nky=int(ky*nk2)
        nkz=int(kz*nk3)
        
        xk=nkx
        yk=nky
        zk=nkz
        
        finek3=nkx*nk2*nk3+nky*nk3+nkz
        for ibnd in range(Libnd):
            #delta=np.real(dkk[i,ibnd,ibnd])/cmax
            EEE=lambdaFS[kindex[1][i],9]
            for ii in range(nbwan):
                index=ii*N*N*N+xk*N*N+yk*N+zk
                #print(kx*nk1,nkx)
                #print(EEE,Ekkxyz[ii,nkx,nky,nkz])
                if abs(Ekkxyz[ii,nkx,nky,nkz]-EEE)<0.01:
                    if nkx==1 and nky==2 and nkz!=0:
                        print(nkx,nky,nkz)
                        print(dkkxyz[ii,nkx,nky,nkz],dkkxyz[ii,nk1-nkx,nk2-nky,nk3-nkz])
                        a=0
                    break
        
        
    finek3=0
    ik=finek3
    
   

            #Ekkxyz[ibnd,nkx,nky,nkz]=EEE

    #
    #frmsf文件中的nk=bxsf文件中的nk-1
    #对frmsf文件中的nk，x=i/(nk)
    #frmsf文件中的nk3对应epw中的nkf3
    effnband=0
    for ii in range(nbwan):
        if np.sum(abs(dkkxyz[ii]))>0.00001:
            effnband=effnband+1
    bxsf=open(bxsffile)
    lead=bxsf.readlines()
    f=open(prefix+str(ne)+'.frmsf','w')
    ii=0
    f.write(str(nk1)+' '+str(nk2)+' '+str(nk3))
    f.write('\n'+'1'+'\n')
    f.write(str(effnband)+'\n')
    for line in lead:
        if ii>=17 and ii<20:
            #print(line)
            f.write(line+'\n')
        ii=ii+1
    for ii in range(int(nbwan)):
        if np.sum(abs(dkkxyz[ii]))>0.00001:
            for xk in range(NN):
                for yk in range(NN):
                    for zk in range(NN):
                        index=ibnd*N*N*N+xk*N*N+yk*N+zk
                        f.write(str(Ekkxyz[ii,xk,yk,zk])+'\n')
                
                
                
    aaaaaaaaa=0
    cmax=abs(np.max(dkkxyz))
    for ii in range(int(nbwan)):
        if np.sum(abs(dkkxyz[ii]))>0.00001:
            for xk in range(NN):
                for yk in range(NN):
                    for zk in range(NN):
                        index=ibnd*N*N*N+xk*N*N+yk*N+zk
                        if abs(Ekkxyz[ibnd,xk,yk,zk])<=1:
                                d=1
                        if xk>0 and xk<NN-1 and yk>0 and yk<NN-1 and zk>0 and zk<NN-1:
                            aaaaaaaaa=dkkxyz[ii,xk,yk,zk]+dkkxyz[ii,xk,yk,zk+1]+dkkxyz[ii,xk,yk,zk-1]+dkkxyz[ii,xk,yk+1,zk]+dkkxyz[ii,xk,yk-1,zk]+dkkxyz[ii,xk+1,yk,zk]+dkkxyz[ii,xk-1,yk,zk]
                        f.write(str(dkkxyz[ii,xk,yk,zk]/cmax)+'\n')    
                        #f.write(str(aaaaaaaaa/cmax)+'\n')
                        #做平滑化之前仍有随机的符号，尚不理解
                        #在自旋方向一致时。对角化得到的delta在一些邻近的k点就出现了符号的突变
                        #怀疑是自旋计算不准
                        #s1=spinorl[ii,xk,yk,zk]
                        #s3=spinorl[ii,xk,yk,zk]
                        #if zk>0:
                        #    ss1=spinorl[ii,xk,yk,zk]
                        #    ss3=spinorl[ii,xk,yk,zk-1]
                        #    r=np.dot(ss1,ss3)
                        #    #print(r)
                        #    if r<0 and r>0.1:
                        #        spinor[ii,xk,yk,zk]=spinor[ii,xk,yk,zk]*-1
                        #f.write(str(spinorll[ii,xk,yk,zk][0])+' '+str(spinorll[ii,xk,yk,zk][1])+' '+str(spinorll[ii,xk,yk,zk][2])+'\n')    
                        
                        #else:
                        #    f.write(str(0)+'\n')
    f.close()
    for xk in range(NN):
        for yk in range(NN):
            for zk in range(NN):
                #print(dkkxyz[2,xk,yk,zk],dkkxyz[2,nk1-xk-1,nk2-yk-1,nk3-zk-1])
                x=0
    #print(dkkxyz)
    #print(effnband)


plotdddd('p')