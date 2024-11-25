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

#bxsf=np.loadtxt('wannier90.bxsf')



wfc_type=1







class irrep_delta():
    def __init__(self, prefix, nk1, nk2, nk3, Ef, ne):
        #parse the paramter of irrep_delta
        self.nk1=nk1
        self.nk2=nk2
        self.nk3=nk3
        self.prefix=prefix
        self.Ef=Ef
        self.ne=ne
        
        #parse the input file
        self.bxsffile=self.prefix+'.bxsf'
        self.Vkkfile= self.prefix+'.lambda_kkq'
        self.datafile= self.prefix+'.lambda_FS'
        
        #get the input data size
        self.Vkkdata=np.loadtxt(self.datafile)
        self.Lk=len(self.Vkkdata[:,3])
        self.Libnd=len(np.unique(self.Vkkdata[:,8]))

        #get k points list
        self.lambdaFS=np.loadtxt(self.datafile)
        self.kindex=np.unique(self.lambdaFS[:,3],return_index=True)
        self.kpoints = self.lambdaFS[self.kindex[1],0:3]
            
        #extract the data we need from input file
        self.nbwan, self.NN, self.Ek, self.Ekkxyz = self.GetElStructure(self.bxsffile, self.nk1+1, self.Ef)
        self.wfck=self.GetWfc(self.Lk, self.nbwan)
        self.dkk0 = self.GetDelta(self.ne, self.Lk, self.Libnd)
        
    #get electron structure information
    def GetElStructure(self, bxsffile, N, Ef):
        bxsf=open(bxsffile)
        lines = bxsf.readlines()
        sum = 0
        Ek = []
        nl=0
        for line in lines:
            nl=nl+1
            if nl >20:
                if line[1]!='B' and line[1]!='E' and line[2]!='E':
                    line=float(line)
                    Ek.append(line)
        
        #calculate thenumber of k_points, band number, minus of electron energy and Ef
        NN=N-1
        nbwan=int(len(Ek)/N/N/N)
        Ekkxyz=np.zeros((nbwan,NN,NN,NN))
        nbwan=int(len(Ek)/N/N/N)
        for ibnd in range(int(nbwan)):
            for xk in range(NN):
                for yk in range(NN):
                    for zk in range(NN):
                        index=ibnd*N*N*N+xk*N*N+yk*N+zk
                        Ekkxyz[ibnd,xk,yk,zk]=Ek[index]-Ef
        return nbwan, NN, Ek, Ekkxyz
    
    #get delta function in /delta directory
    def GetDelta(self, ne, Lk, Libnd):
        eevvv = np.loadtxt('delta'+'/'+'dkk'+str(ne)+'.txt').view(complex).reshape(-1)
        lll=abs(np.sqrt((eevvv**2).sum(axis=0)))
        dkk0=np.zeros((Lk,Libnd,Libnd),dtype=np.complex64)
        for i in range(Lk):
            for iii in range(int(Libnd/2)):
                for iibnd in range(2):
                    for jjbnd in range(2):
                        dkk0[i][2*iii+iibnd, 2*iii+jjbnd] = eevvv[i * Libnd * Libnd +iii*4+ iibnd * 2 + jjbnd]
        return dkk0
    
    #get electron wavefunction in wannier basis in wfc333.dat
    def GetWfc(self, Lk, nbwan):
        d1=np.loadtxt('wfc.dat')
        wfck=np.zeros((nbwan,Lk,nbwan),dtype=np.complex64)
        LLwfc=len(d1)
        for i in range(LLwfc):
            ik=int(d1[i,0]-1)
            ibnd=int(d1[i,1]-1)
            jbnd=int(d1[i,2]-1)
            imode=int(d1[i,3]-1)
            if imode==0:
                wfck[ibnd,ik,jbnd]=d1[i,4]+1j*d1[i,5]
        return wfck
        
    #align delta function to electron structure
    def AlignDeltaToElStructure(self, Lk, nk1, nk2, nk3, Libnd, nbwan, NN, Ekkxyz, wfck, dkk0, lambdaFS, kindex):
        dkkxyz=np.zeros((nbwan,NN,NN,NN))
        spin_xyz = np.zeros((nbwan,NN,NN,NN,3))
        cmax=max(abs(dkk0[:,0,0]))
        for i in range(Lk):
            #get k points index and position in lambdaFS(include all delta function k points we calculated)
            ik=kindex[1][i]
            kk1=int(lambdaFS[ik,3])-1
            kkk1=int(lambdaFS[ik,4])-1
            
            kx=lambdaFS[kindex[1][i],0]
            ky=lambdaFS[kindex[1][i],1]
            kz=lambdaFS[kindex[1][i],2]
            
            nkx=int(kx*nk1)
            nky=int(ky*nk2)
            nkz=int(kz*nk3)

            sum=0
            for ii in range(int(Libnd/2)):
                
                EEE=lambdaFS[kindex[1][i],9]
                for iii in range(nbwan):
                    #align electruture to delta function
                    if abs(Ekkxyz[iii,nkx,nky,nkz]-EEE)<0.002 and iii%2==0:
                        #get -k position
                        nkx_inverse=nkx
                        nky_inverse=nky
                        nkz_inverse=nkz
                        if nkz!=0:
                            nkz_inverse=nk1-nkz
                        if nky!=0:
                            nky_inverse=nk2-nky
                        if nkx!=0:
                            nkx_inverse=nk1-nkx

                        #diagnolized the delta function
                        diakk, diag_delta_k = self.DiagDelta(dkk0[kk1], ii)
                        diakkk, diag_delta_inverse_k = self.DiagDelta(dkk0[kkk1], ii)

                        #transform basis to diagnolized delta function
                        wfckk=[wfck[iii,kk1],wfck[iii+1,kk1]]
                        wfckk0, wfckk1 = self.TranWfc(wfckk, diakk)
                        
                        wfckkk=[wfck[iii,kkk1],wfck[iii+1,kkk1]]
                        wfckkk0, wfckkk1 = self.TranWfc(wfckkk, diakkk)
                        
                        #save spin data
                        [sx, sy, sz] = calc_spin_wfc(wfckk0)
                        [sx1, sy1, sz1] = calc_spin_wfc(wfckkk0)
                        spin_xyz[iii,nkx,nky,nkz] = [sx, sy, sz]
                        if iii%2==0:
                            dkkxyz[iii,nkx,nky,nkz]=np.real(diag_delta_k[0,0])
                            if sx<0:
                                dkkxyz[iii,nkx,nky,nkz]=-np.real(diag_delta_k[0,0])
        return dkkxyz, spin_xyz
    
    #diagnolized the delta function
    def DiagDelta(self, dkk, ii):
        delta_k = np.array([[dkk[2*ii, 2*ii], dkk[2*ii, 2*ii+1]], [dkk[2*ii+1, 2*ii], dkk[2*ii+1, 2*ii+1]]])
        delta_k =delta_k+np.conjugate(np.transpose(delta_k))
        [eigkk, diakk] = np.linalg.eig(delta_k)
        diag_delta_k =  np.dot(np.dot(np.conjugate(np.transpose(diakk)), delta_k), diakk)
        return diakk, diag_delta_k
    #transform the wavefunction of electron to the basis of diaglized delta function
    def TranWfc(self, wfckk, diakk):
        
        wfckk_diag=np.dot(np.conjugate(np.transpose(diakk)),wfckk)
        wfckk0=wfckk_diag[0,:]
        wfckk1=wfckk_diag[1,:]
        return wfckk0, wfckk1
    
    #delta fuction of /delta only finite when abs(Ek-Ef)<fthick, we need interpolate delta function to full k space for plot the delta function
    def InterpolateDelta(self, nk1, nk2, nk3, nbwan, NN, dkkxyz):
        #initialize interpolated delta function
        dkkxyz_interpolate = dkkxyz
        
        #interpolate delta function of every band
        for i in range(nbwan):
            rawdata=[]
            #points=np.array([[0,0,0]])
            x=[]
            y=[]
            z=[]
            

            #generate input data for interpolate function Rbf
            for nkx in range(NN):
                for nky in range(NN):
                    for nkz in range(NN):

                        if dkkxyz[i,nkx,nky,nkz]>0.00001 or dkkxyz[i,nkx,nky,nkz]<-0.00001:
                            rawdata.append(dkkxyz[i,nkx,nky,nkz])
                            x.append(nkx)
                            y.append(nky)
                            z.append(nkz)
            points=np.array(np.transpose([x,y,z]))
            #interpolate delta function to dkkxyz_interpolate
            if len(points)>=1:
                fun = Rbf(x,y,z,rawdata,method='linear',smooth=0)
                for nkx in range(nk1):
                    for nky in range(nk2):
                        for nkz in range(nk3):
                            dkkxyz_interpolate[i,nkx,nky,nkz]=fun(nkx,nky,nkz)
        return dkkxyz_interpolate
    
    #dump delta function in .frmsf format for plot figure
    def SaveToFrmsf(self, bxsffile, nk1, nk2, nk3, NN, nbwan, prefix, ne, Ekkxyz, dkkxyz):
        #write the stucture information, k points information, band information to .frmsf file by .bxsf file
        effnband=0
        for ii in range(nbwan):
            print(ii)
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
            
        #calculate electron structure information to .frmsf file
        for ii in range(int(nbwan)):
            print(ii)
            if np.sum(abs(dkkxyz[ii]))>0.00001:
                for xk in range(NN):
                    for yk in range(NN):
                        for zk in range(NN):
                            f.write(str(Ekkxyz[ii,xk,yk,zk])+'\n')
                    
        #calculate delta function as color information to .frmsf file
        cmax=abs(np.max(dkkxyz))
        for ii in range(int(nbwan)):
            if np.sum(abs(dkkxyz[ii]))>0.00001:
                for xk in range(NN):
                    for yk in range(NN):
                        for zk in range(NN):
                            f.write(str(dkkxyz[ii,xk,yk,zk]/cmax)+'\n')    
        f.close()
    def CheckRepC2x(self, dkkxyz, nbwan, NN, Lk, kindex, lambdaFS, Ekkxyz):
        trans_matrix = [[0, -1, 1], [0, -1, 0], [1, -1, 0]]
        #trans_matrix = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
        for ibnd in range(nbwan):
            for i in range(Lk):
                ik=kindex[1][i]
                kk1=int(lambdaFS[ik,3])-1
                kkk1=int(lambdaFS[ik,4])-1
                
                kx=lambdaFS[kindex[1][i],0]
                ky=lambdaFS[kindex[1][i],1]
                kz=lambdaFS[kindex[1][i],2]
                
                nkx=int(kx*NN+10**-6)
                nky=int(ky*NN+10**-6)
                nkz=int(kz*NN+10**-6)
                k_vector = [kx, ky, kz]
                if abs(dkkxyz[ibnd, nkx, nky, nkz]) > 0.000001 and ibnd%2==0: 
                        
                    nkx_trans, nky_trans, nkz_trans = self.TransKvec(trans_matrix, k_vector, NN)
                    print("n",nkx, nky, nkz)
                    print(nkx_trans, nky_trans, nkz_trans)
                    #print(int(nkx_trans+10**-6), int(nky_trans+10**-6), int(nkz_trans+10**-6))
                    print(Ekkxyz[ibnd, nkx, nky, nkz], Ekkxyz[ibnd, nkx_trans, nky_trans, nkz_trans])
                    print(dkkxyz[ibnd, nkx, nky, nkz], dkkxyz[ibnd, nkx_trans, nky_trans, nkz_trans])
                    print(self.spin_xyz[ibnd, nkx, nky, nkz], self.spin_xyz[ibnd, nkx_trans, nky_trans, nkz_trans])
    def CheckSolution(self, dkkxyz, nbwan, NN, Lk, kindex, lambdaFS, Ekkxyz):
        trans_matrix = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]
        #trans_matrix = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
        for ibnd in range(nbwan):
            for i in range(Lk):
                ik=kindex[1][i]
                kk1=int(lambdaFS[ik,3])-1
                kkk1=int(lambdaFS[ik,4])-1
                
                kx=lambdaFS[kindex[1][i],0]
                ky=lambdaFS[kindex[1][i],1]
                kz=lambdaFS[kindex[1][i],2]
                
                nkx=int(kx*NN+10**-6)
                nky=int(ky*NN+10**-6)
                nkz=int(kz*NN+10**-6)
                k_vector = [kx, ky, kz]
                if abs(dkkxyz[ibnd, nkx, nky, nkz]) > 0.000001 and ibnd%2==0: 
                        
                    nkx_trans, nky_trans, nkz_trans = self.TransKvec(trans_matrix, k_vector, NN)
                    print("nk:",nkx, nky, nkz)
                    print("nk_invers:",nkx_trans, nky_trans, nkz_trans)
                    print("delta(n,k) = ", dkkxyz[ibnd, nkx, nky, nkz], "delta(n,-k) = ",dkkxyz[ibnd, nkx_trans, nky_trans, nkz_trans])
                    print("spin(n,k) = ", self.spin_xyz[ibnd, nkx, nky, nkz], "spin(n,-k) = ", self.spin_xyz[ibnd, nkx_trans, nky_trans, nkz_trans])
    def TransKvec(self, trans_matrix, k_vector, NN):
        print("k",k_vector)
        for i in range(len(k_vector)):
        
            if k_vector[i] >= 0.5:
                k_vector[i] = k_vector[i] - 1
        k_trans_vector = np.dot(trans_matrix, k_vector)
        for i in range(len(k_vector)):
            if k_trans_vector[i] < 0:
                k_trans_vector[i] = k_trans_vector[i] + 1
        
        print(k_trans_vector)
        return int(k_trans_vector[0]*NN+10**-6), int(k_trans_vector[1]*NN+10**-6), int(k_trans_vector[2]*NN+10**-6)
        
    def RunAll(self):
        
        #align delta function to electron structure
        dkkxyz, self.spin_xyz = self.AlignDeltaToElStructure(self.Lk, self.nk1, self.nk2, self.nk3, self.Libnd, self.nbwan, self.NN,  self.Ekkxyz, self.wfck, self.dkk0, self.lambdaFS, self.kindex)
        
        #check whether the gap function is physical meaningful
        self.CheckSolution(dkkxyz, self.nbwan, self.NN, self.Lk, self.kindex, self.lambdaFS, self.Ekkxyz)
        
        #check what the character of C2x operation in this gap function
        self.CheckRepC2x(dkkxyz, self.nbwan, self.NN, self.Lk, self.kindex, self.lambdaFS, self.Ekkxyz)
        
        #delta fuction of /delta only finite when abs(Ek-Ef)<fthick, we need interpolate delta function to full k space for plot the delta function
        dkkxyz = self.InterpolateDelta(self.nk1, self.nk2, self.nk3, self.nbwan, self.NN, dkkxyz)
    
        
        #dump delta function in .frmsf format for plot figure
        self.SaveToFrmsf(self.bxsffile, self.nk1, self.nk2, self.nk3, self.NN, self.nbwan, self.prefix, self.ne, self.Ekkxyz, dkkxyz)
            #Ekkxyz[ibnd,nkx,nky,nkz]=EEE



irrep = irrep_delta("Pb", 20, 20, 20, 11.2876, 4)
irrep.RunAll()