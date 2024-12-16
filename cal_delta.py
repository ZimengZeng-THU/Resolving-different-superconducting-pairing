#!~/.conda/envs/nlp/bin/python3.7
import numpy as np
import matplotlib.pyplot as plt
import math as mmm
import numpy.linalg as la
from scipy.sparse.linalg.eigen.arpack import eigs as largest_eigsh
import csv
import os

plt.rc("font", family='Times New Roman')


class cal_delta():
    def __init__(self, prefix, nk1, nk2, nk3, nmode):
        #parse class parameters
        self.NE=36
        self.nk1=nk1
        self.nk2=nk2
        self.nk3=nk3
        self.nmode=nmode
        
        #get input data and create output format
        self.dosef0, self.dosef1, self.dosef_dos, self.Vkkdata1, self.Vkk_s, self.Vkk_t, self.dkk_s, self.dkk_t = self.GetInput(prefix, nmode)
    
    #get the size of input data(number of kpoints, number of bands, size of .kkq file)
    def GetInputSize(self, Vkkdata, lambdaFS):
        LVkk = len(Vkkdata)
        kkkk = np.unique(Vkkdata[:, 0])
        aaaa = np.unique(Vkkdata[:, 4])
        kindex = np.unique(lambdaFS[:, 3], return_index=True)
        Lk = len(kindex[0])
        Libnd = len(aaaa)
        Lkk = np.zeros((Lk))
        return LVkk, kkkk, Lk, Libnd, Lkk
    
    #get all data we need by input file
    def GetInput(self, prefix, nmode):
        Vkkfile = prefix+'.lambda_kkq'
        Vkkdata = np.loadtxt(Vkkfile)
        datafile = prefix+'.lambda_FS'
        lambdaFS = np.loadtxt(datafile)
        LVkk, kkkk, Lk, Libnd, Lkk = self.GetInputSize(Vkkdata, lambdaFS)
        self.lambdaFS = lambdaFS
        self.LVkk = LVkk
        self.Lk = Lk
        self.Lkk = Lkk
        self.Libnd = Libnd
        
        #Generate containers that store intermediate variables
        Vkkdata1 = np.zeros((Lk, Lk, Libnd, Libnd, nmode), dtype=np.complex64)
        dosef0 = np.zeros((Lk, Libnd))
        dosef_dos = np.zeros((Lk, Libnd))
        dosef1 = np.zeros((Lk, Libnd))
        Vkk_s = np.zeros((Lk * int(Libnd/2) * int(Libnd/2), Lk * int(Libnd/2) * int(Libnd/2)), dtype=np.complex64)
        Vkk_t = np.zeros((Lk * int(Libnd/2) * int(Libnd/2) * 3, Lk * int(Libnd/2) * int(Libnd/2) * 3), dtype=np.complex64)
        dkk_s = np.zeros((Lk, int(Libnd/2), int(Libnd/2)), dtype=np.complex64)
        dkk_t = np.zeros((Lk, int(Libnd/2), int(Libnd/2)), dtype=np.complex64)
        
        ii=0
        sum1=0
        #get dos data, el-ph matrix data
        while ii < LVkk:
            ik = int(Vkkdata[ii, 0]) - 1
            ikq = int(Vkkdata[ii, 2]) - 1
            sk1 = int(Vkkdata[ii, 4]) - 1
            skq1 = int(Vkkdata[ii, 5]) - 1
            iskq10 = skq1 - 2 * int(skq1 / 2)
            iskq11 = 2 * int(skq1 / 2)
            imode = int(Vkkdata[ii, 8]) - 1
        
            dosef0[ik, sk1] = 1 / self.nk1 / self.nk2 / self.nk3 * Vkkdata[ii, 10]
            dosikq = dosef0[ikq, iskq11 + iskq10] / 2 + dosef0[ikq, iskq11 + 1 - iskq10] / 2
            dosef1[ikq, skq1] = dosikq
            dosef_dos[ik,sk1]=1/self.nk1/self.nk2/self.nk3*Vkkdata[ii,12]
        
            g1 = Vkkdata[ii, 6] + 1j * Vkkdata[ii, 7]
        
            if Vkkdata[ii, 9] > 0.0000001:
                Vkkdata1[ik, ikq, sk1, skq1, imode] = g1 / np.sqrt(Vkkdata[ii, 9])
                sum1 = sum1 + abs(g1 * g1 / (self.nk1 * self.nk2 * self.nk3) / (self.nk1 * self.nk2 * self.nk3) * Vkkdata[ii, 10] * Vkkdata[ii, 11]) / \
                   Vkkdata[ii, 9]
            ii=ii+1
        print(sum1, sum1 / np.sum(dosef_dos) * 2)
        return dosef0, dosef1, dosef_dos, Vkkdata1, Vkk_s, Vkk_t, dkk_s, dkk_t
        
    #calculate interaction matrix
    def CalVkk(self, Lk, Libnd, Vkkdata1, dosef1, nmode):
        Vkk = np.zeros((Lk * Libnd * Libnd, Lk * Libnd * Libnd), dtype=np.complex64)
        Vkk_s = np.zeros((Lk * int(Libnd/2) * int(Libnd/2), Lk * int(Libnd/2) * int(Libnd/2)), dtype=np.complex64)
        #Vkk_s = np.zeros((Lk * Libnd * Libnd, Lk * Libnd * Libnd), dtype=np.complex64)
        Vkk_t = np.zeros((Lk * int(Libnd/2) * int(Libnd/2) * 3, Lk * int(Libnd/2) * int(Libnd/2) * 3), dtype=np.complex64)
        trans_matrix = [[1/np.sqrt(2),0,-1/np.sqrt(2),0],[0,1,0,0],[1/np.sqrt(2),0,1/np.sqrt(2),0],[0,0,0,1]]
        for ik in range(Lk):
            for ikq in range(Lk):
                for sk1 in range(Libnd):
                    for sk2 in range(Libnd):
                        for skq1 in range(Libnd):
                            for skq2 in range(Libnd):
                                for imode in range(nmode):
                                    spin_factor_k = 1
                                    if sk2%2==0:
                                        spin_factor_k = 1
                                    else:
                                        spin_factor_k = -1
                                    spin_factor_kq = 1
                                    if skq2%2==0:
                                        spin_factor_kq = 1
                                    else:
                                        spin_factor_kq = -1
                                    Vkk[ik*Libnd*Libnd+sk1*Libnd+sk2,ikq*Libnd*Libnd+(skq1)*Libnd+skq2] = \
                                            Vkk[ik*Libnd*Libnd+sk1*Libnd+sk2,ikq*Libnd*Libnd+(skq1)*Libnd+skq2] \
                                            - ( Vkkdata1[ik,ikq,sk1,skq1,imode]*np.conjugate(Vkkdata1[ik,ikq,sk2,skq2,imode])) \
                                            *dosef1[ikq,skq1]*2
                        
        if 1 == 1:
            for ik in range(Lk):
                for ikq in range(Lk):
                    V00=Vkk[ik*Libnd*Libnd+0*Libnd+0,ikq*Libnd*Libnd+(0)*Libnd+0]
                    V01=Vkk[ik*Libnd*Libnd+0*Libnd+0,ikq*Libnd*Libnd+(0)*Libnd+1]
                    V02=Vkk[ik*Libnd*Libnd+0*Libnd+0,ikq*Libnd*Libnd+(1)*Libnd+1]
                    V03=Vkk[ik*Libnd*Libnd+0*Libnd+0,ikq*Libnd*Libnd+(1)*Libnd+0]
                    
                    V10=Vkk[ik*Libnd*Libnd+0*Libnd+1,ikq*Libnd*Libnd+(0)*Libnd+0]
                    V11=Vkk[ik*Libnd*Libnd+0*Libnd+1,ikq*Libnd*Libnd+(0)*Libnd+1]
                    V12=Vkk[ik*Libnd*Libnd+0*Libnd+1,ikq*Libnd*Libnd+(1)*Libnd+1]
                    V13=Vkk[ik*Libnd*Libnd+0*Libnd+1,ikq*Libnd*Libnd+(1)*Libnd+0]
                    
                    V30=Vkk[ik*Libnd*Libnd+1*Libnd+0,ikq*Libnd*Libnd+(0)*Libnd+0]
                    V31=Vkk[ik*Libnd*Libnd+1*Libnd+0,ikq*Libnd*Libnd+(0)*Libnd+1]
                    V32=Vkk[ik*Libnd*Libnd+1*Libnd+0,ikq*Libnd*Libnd+(1)*Libnd+1]
                    V33=Vkk[ik*Libnd*Libnd+1*Libnd+0,ikq*Libnd*Libnd+(1)*Libnd+0]
                    
                    V20=Vkk[ik*Libnd*Libnd+1*Libnd+1,ikq*Libnd*Libnd+(0)*Libnd+0]
                    V21=Vkk[ik*Libnd*Libnd+1*Libnd+1,ikq*Libnd*Libnd+(0)*Libnd+1]
                    V22=Vkk[ik*Libnd*Libnd+1*Libnd+1,ikq*Libnd*Libnd+(1)*Libnd+1]
                    V23=Vkk[ik*Libnd*Libnd+1*Libnd+1,ikq*Libnd*Libnd+(1)*Libnd+0]
                    Vkk4=np.array([[V00,V01,V02,V03],[V10,V11,V12,V13],[V20,V21,V22,V23],[V30,V31,V32,V33]])
                    #print("before",Vkk4[0,1],Vkk4[1,0])
                    Vkk4=np.dot(np.dot(trans_matrix,Vkk4), np.transpose(trans_matrix))
                    #print("after",Vkk4[0,1],Vkk4[1,0])
                    #print('\n')
                    Vkk_s[ik,ikq] = Vkk4[0,0]
                                         
        if 1 == 1:
            for ik in range(Lk):
                for ikq in range(Lk):
                    V00=Vkk[ik*Libnd*Libnd+0*Libnd+0,ikq*Libnd*Libnd+(0)*Libnd+0]
                    V01=Vkk[ik*Libnd*Libnd+0*Libnd+0,ikq*Libnd*Libnd+(0)*Libnd+1]
                    V02=Vkk[ik*Libnd*Libnd+0*Libnd+0,ikq*Libnd*Libnd+(1)*Libnd+1]
                    V03=Vkk[ik*Libnd*Libnd+0*Libnd+0,ikq*Libnd*Libnd+(1)*Libnd+0]
                    
                    V10=Vkk[ik*Libnd*Libnd+0*Libnd+1,ikq*Libnd*Libnd+(0)*Libnd+0]
                    V11=Vkk[ik*Libnd*Libnd+0*Libnd+1,ikq*Libnd*Libnd+(0)*Libnd+1]
                    V12=Vkk[ik*Libnd*Libnd+0*Libnd+1,ikq*Libnd*Libnd+(1)*Libnd+1]
                    V13=Vkk[ik*Libnd*Libnd+0*Libnd+1,ikq*Libnd*Libnd+(1)*Libnd+0]
                    
                    V30=Vkk[ik*Libnd*Libnd+1*Libnd+0,ikq*Libnd*Libnd+(0)*Libnd+0]
                    V31=Vkk[ik*Libnd*Libnd+1*Libnd+0,ikq*Libnd*Libnd+(0)*Libnd+1]
                    V32=Vkk[ik*Libnd*Libnd+1*Libnd+0,ikq*Libnd*Libnd+(1)*Libnd+1]
                    V33=Vkk[ik*Libnd*Libnd+1*Libnd+0,ikq*Libnd*Libnd+(1)*Libnd+0]
                    
                    V20=Vkk[ik*Libnd*Libnd+1*Libnd+1,ikq*Libnd*Libnd+(0)*Libnd+0]
                    V21=Vkk[ik*Libnd*Libnd+1*Libnd+1,ikq*Libnd*Libnd+(0)*Libnd+1]
                    V22=Vkk[ik*Libnd*Libnd+1*Libnd+1,ikq*Libnd*Libnd+(1)*Libnd+1]
                    V23=Vkk[ik*Libnd*Libnd+1*Libnd+1,ikq*Libnd*Libnd+(1)*Libnd+0]
                    Vkk4=np.array([[V00,V01,V02,V03],[V10,V11,V12,V13],[V20,V21,V22,V23],[V30,V31,V32,V33]])
                    #print("before",Vkk4[0,1],Vkk4[1,0])
                    Vkk4=np.dot(np.dot(trans_matrix,Vkk4), np.transpose(trans_matrix))
                    #print("after",Vkk4[0,1],Vkk4[1,0])
                    #print('\n')

                    for imode in range(nmode):
                        spin_factor_k = 1
                        if sk2%2==0:
                            spin_factor_k = 1
                        else:
                            spin_factor_k = -1
                        spin_factor_kq = 1
                        if skq2%2==0:
                            spin_factor_kq = 1
                        else:
                            spin_factor_kq = -1
                        i=0
                        j=0
                        Vkk_t[ik*3+i,ikq*3+j] = Vkk_t[ik*3+i,ikq*3+j] \
                                - ( Vkkdata1[ik,ikq,0,0,imode]*np.conjugate(Vkkdata1[ik,ikq,1,1,imode])) \
                                *dosef1[ikq,skq1]*2


                        i=0
                        j=1
                        Vkk_t[ik*3+i,ikq*3+j] = Vkk_t[ik*3+i,ikq*3+j] \
                                - ( -Vkkdata1[ik,ikq,0,0,imode]*np.conjugate(Vkkdata1[ik,ikq,1,0,imode])
                                   +Vkkdata1[ik,ikq,0,1,imode]*np.conjugate(Vkkdata1[ik,ikq,1,1,imode])
                                   ) \
                                *dosef1[ikq,skq1]*2 / np.sqrt(2)

                        i=0
                        j=2
                        Vkk_t[ik*3+i,ikq*3+j] = Vkk_t[ik*3+i,ikq*3+j] \
                                - ( -Vkkdata1[ik,ikq,0,1,imode]*np.conjugate(Vkkdata1[ik,ikq,1,0,imode])) \
                                *dosef1[ikq,skq1]*2
                        
                        
                        i=1
                        j=0
                        Vkk_t[ik*3+i,ikq*3+j] = Vkk_t[ik*3+i,ikq*3+j] \
                                - ( -Vkkdata1[ik,ikq,0,0,imode]*np.conjugate(Vkkdata1[ik,ikq,0,1,imode]) \
                                    +Vkkdata1[ik,ikq,1,0,imode]*np.conjugate(Vkkdata1[ik,ikq,1,1,imode])
                                   ) \
                                *dosef1[ikq,skq1]*2 / np.sqrt(2)
                        
                        i=1
                        j=1
                        Vkk_t[ik*3+i,ikq*3+j] = Vkk_t[ik*3+i,ikq*3+j] \
                                - ( Vkkdata1[ik,ikq,0,0,imode]*np.conjugate(Vkkdata1[ik,ikq,0,0,imode])
                                   -Vkkdata1[ik,ikq,0,1,imode]*np.conjugate(Vkkdata1[ik,ikq,0,1,imode])
                                   -Vkkdata1[ik,ikq,1,0,imode]*np.conjugate(Vkkdata1[ik,ikq,1,0,imode])
                                   +Vkkdata1[ik,ikq,1,1,imode]*np.conjugate(Vkkdata1[ik,ikq,1,1,imode])
                                   ) \
                                *dosef1[ikq,skq1]*2 / 2
                                
                        i=1
                        j=2
                        Vkk_t[ik*3+i,ikq*3+j] = Vkk_t[ik*3+i,ikq*3+j] \
                                - ( Vkkdata1[ik,ikq,0,1,imode]*np.conjugate(Vkkdata1[ik,ikq,0,0,imode])
                                   -Vkkdata1[ik,ikq,1,1,imode]*np.conjugate(Vkkdata1[ik,ikq,1,0,imode])
                                
                                   ) \
                                *dosef1[ikq,skq1]*2 / np.sqrt(2)
                                
                                
                        i=2
                        j=0
                        Vkk_t[ik*3+i,ikq*3+j] = Vkk_t[ik*3+i,ikq*3+j] \
                                - ( -Vkkdata1[ik,ikq,1,0,imode]*np.conjugate(Vkkdata1[ik,ikq,0,1,imode])
                                   ) \
                                *dosef1[ikq,skq1]*2
                        
                        
                        i=2
                        j=1
                        Vkk_t[ik*3+i,ikq*3+j] = Vkk_t[ik*3+i,ikq*3+j] \
                               - ( -Vkkdata1[ik,ikq,1,1,imode]*np.conjugate(Vkkdata1[ik,ikq,0,1,imode])
                                  +Vkkdata1[ik,ikq,1,0,imode]*np.conjugate(Vkkdata1[ik,ikq,0,0,imode])
                               
                                  ) \
                               *dosef1[ikq,skq1]*2 / np.sqrt(2)
                        
                        i=2
                        j=2
                        Vkk_t[ik*3+i,ikq*3+j] = Vkk_t[ik*3+i,ikq*3+j] \
                                - ( Vkkdata1[ik,ikq,1,1,imode]*np.conjugate(Vkkdata1[ik,ikq,0,0,imode])
                                   ) \
                                *dosef1[ikq,skq1]*2
                        
                                            
        return Vkk,Vkk_s,Vkk_t
    
    #calculate eigenvector of Vkk
    def CalEigen(self, Vkk, NE):
        eaa, evv = largest_eigsh(Vkk, NE, which='SR')
        evals, evecs = largest_eigsh(Vkk, NE, which='SR')
        return evals,evecs
    
    #the eigenvector of Vkk carry a random phase, this function get this function
    def GetPhase(self, dkk, Lk):
        for i in range(Lk):
            if abs(dkk[i][0, 0]) > 0:
                a = dkk[i][0, 0] / abs(dkk[i][0, 0])
                #print(abs(dkk[i][0, 0]))
                break
        return a
    
    #save all data we need
    def SaveAll(self,dkk_s, dkk_t, dkk):
        for vv in range(self.NE):
            self.WriteDelta(dkk_s[vv], "s", vv, self.Lk, self.Libnd)
            self.WriteDelta(dkk_t[vv], "t", vv, self.Lk, self.Libnd)
            self.WriteDelta(dkk[vv], "", vv, self.Lk, self.Libnd)
    #dump delta function data
    def WriteDelta(self, dkk, parity, aaaaa, Lk, Libnd):
        evv = np.zeros((Lk *Libnd *Libnd), dtype = np.complex64)
        #create delta directory
        if not os.path.exists("delta"):
            os.mkdir("delta")
        #save delta function to vector
        for i in range(Lk):
            deltakk = np.array([[dkk[i][0, 0], -dkk[i][0, 1]], [-dkk[i][1, 0], dkk[i][1, 1]]])
            for ibnd in range(2):
                for jbnd in range(2):
                    evv[i * Libnd * Libnd +ibnd*2 + jbnd] = dkk[i][ibnd, jbnd]
        #save vector
        np.savetxt('delta/' + 'dkk' + parity + str(aaaaa) + '.txt', np.column_stack([evv.real, evv.imag]))
        
    #Change the storage format of eigenvectors to store the four components
    #(singlet:\sigma_0; triplet:\sigma_1, \sigma_2, \sigma_3) of each k-point as a 2 * 2 matrix
    def SaveToDelta(self, NE, evecs, parity, Lk, Libnd):
        dkk = np.zeros((NE, Lk, Libnd, Libnd), dtype=np.complex64)
        for vv in range(NE):
            phase = 1
            for i in range(Lk):
                if parity == "s":
                    dkk[vv][i][0, 0] = evecs[i, vv] / np.sqrt(2)
                    dkk[vv][i][1, 1] = evecs[i, vv] / np.sqrt(2)
                if parity == "t":
                    dkk[vv][i][0, 0] = evecs[i * 3 + 1, vv] / np.sqrt(2)
                    dkk[vv][i][1, 1] = -evecs[i * 3 + 1, vv] / np.sqrt(2)
                    dkk[vv][i][0, 1] = -evecs[i * 3 + 0, vv]
                    dkk[vv][i][1, 0] = evecs[i * 3 + 2, vv]
                if parity == "all":
                    for sk1 in range(2):
                        for skq1 in range(2):
                            dkk[vv][i][sk1, skq1] = evecs[i * Libnd *Libnd +sk1*Libnd +skq1, vv]
                #if i == 0 and vv == 1:
                #    print(vv, self.dkk[vv][i])
            phase = self.GetPhase(dkk[vv], self.Lk)
            for i in range(Lk):
                dkk[vv][i] = dkk[vv][i] / phase
                if i == 0:
                    print(vv, dkk[vv][i])
                    #print(evecs[0,vv],evecs[1,vv],evecs[2,vv])
        return dkk

    #normalize our results by the density of states of EPW calculation
    def NormVal(self, evals, dosef_k_dos, dosef_k_epc):
        dosef_single_spin_dos = sum(dosef_k_dos[:,0]) / 2 + sum(dosef_k_dos[:,1]) / 2
        dosef_single_spin_epc = sum(dosef_k_epc[:,0]) / 2 + sum(dosef_k_epc[:,1]) / 2
        evals = evals * dosef_single_spin_epc / dosef_single_spin_dos
        return evals
    
    def RunAll(self):
        
        #calculate the interaction matrix
        Vkk,Vkk_s,Vkk_t = self.CalVkk(self.Lk, self.Libnd, self.Vkkdata1, self.dosef1, self.nmode)
        
        #setting how many eigenstates we calculate
        NE=36
        
        #calculate the eigenvector of interaction matrix Vkk
        evals_s, evecs_s = self.CalEigen(Vkk_s, NE)
        evals_t, evecs_t = self.CalEigen(Vkk_t, NE)
        evals, evecs = self.CalEigen(Vkk, NE)
        #because the define of density of states of dos calculation and epc calculation is different,
        #we normalize our results by the density of states of EPW calculation
        evals_s = self.NormVal(evals_s, self.dosef_dos, self.dosef1)
        evals_t = self.NormVal(evals_t, self.dosef_dos, self.dosef1)
        
        #save our delta results
        dkk_s = self.SaveToDelta(NE, evecs_s, "s", self.Lk, self.Libnd)
        dkk_t = self.SaveToDelta(NE, evecs_t, "t", self.Lk, self.Libnd)
        dkk = self.SaveToDelta(NE, evecs, "all", self.Lk, self.Libnd)
        self.SaveAll(dkk_s, dkk_t, dkk)

        print(evals_s)
        print(evals_t)
#cal_delta(prefix, nk1,nk2, nk3, nmode)
new = cal_delta("pb", 10, 10, 10, 3)
new.RunAll()
