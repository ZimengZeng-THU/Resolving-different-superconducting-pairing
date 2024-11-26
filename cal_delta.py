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
        self.dosef0, self.dosef1, self.dosef_dos, self.Vkkdata1, self.Vkk, self.dkk = self.GetInput(prefix, nmode)
    
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
        Vkk = np.zeros((Lk * Libnd * Libnd, Lk * Libnd * Libnd), dtype=np.complex64)
        dkk = np.zeros((Lk, Libnd, Libnd), dtype=np.complex64)
        
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
        return dosef0, dosef1, dosef_dos, Vkkdata1, Vkk, dkk
        
    #calculate interaction matrix
    def CalVkk(self, Lk, Libnd, Vkk, Vkkdata1, dosef1, nmode):
        for ik in range(Lk):
            for ikq in range(Lk):
                for sk1 in range(Libnd):
                    for sk2 in range(Libnd):
                        for skq1 in range(Libnd):
                            for skq2 in range(Libnd):
                                for imode in range(nmode):
                                    Vkk[ik*Libnd*Libnd+sk1*Libnd+sk2,ikq*Libnd*Libnd+(skq1)*Libnd+skq2]=Vkk[ik*Libnd*Libnd+sk1*Libnd+sk2,ikq*Libnd*Libnd+(skq1)*Libnd+skq2]-Vkkdata1[ik,ikq,sk1,skq1,imode]*np.conjugate(Vkkdata1[ik,ikq,sk2,skq2,imode])*dosef1[ikq,skq1]*2
                                phase=0
        return Vkk
    
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
    def SaveAll(self):
        for vv in range(self.NE):
            self.WriteDelta(self.dkk[vv], vv, self.Lk, self.Libnd)
            
    #dump delta function data
    def WriteDelta(self, dkk, aaaaa, Lk, Libnd):
        eev = np.zeros((Lk * Libnd * Libnd), dtype=np.complex64)
        #create delta directory
        if not os.path.exists("delta"):
            os.mkdir("delta")
        #save delta function to vector
        for i in range(Lk):
            for ii in range(int(Libnd/2)):
                deltakk = np.array([[dkk[i][2*ii, ii], dkk[i][2*ii, 2*ii+1]], [dkk[i][2*ii+1, 2*ii], dkk[i][2*ii+1, 2*ii+1]]])
                #print(deltakk)
                [eigkk, diakk] = np.linalg.eig(deltakk)

                if np.linalg.det(diakk) == -1:
                    signkk = -1
                else:
                    signkk = 1
                newdeltakk = signkk * np.dot(np.dot(np.conjugate(np.transpose(diakk)), deltakk), diakk)
                for ibnd in range(2):
                    for jbnd in range(2):
                        eev[i * Libnd * Libnd + ii*4+ibnd*2 + jbnd] = dkk[i][ibnd, jbnd]
        #save vector
        np.savetxt('delta/' + 'dkk' + str(aaaaa) + '.txt', np.column_stack([eev.real, eev.imag]))
        
    #Change the storage format of eigenvectors to store the four components
    #(singlet:\sigma_0; triplet:\sigma_1, \sigma_2, \sigma_3) of each k-point as a 2 * 2 matrix
    def SaveToDelta(self, NE, evecs, Lk, Libnd):
        self.dkk = np.zeros((NE, Lk, Libnd, Libnd), dtype=np.complex64)
        for vv in range(NE):
            phase = 1
            for i in range(Lk):
                for iibnd in range(Libnd):
                    for jjbnd in range(Libnd):
                        self.dkk[vv][i][iibnd, jjbnd] = evecs[i * Libnd * Libnd + iibnd * Libnd + jjbnd, vv]
            phase = self.GetPhase(self.dkk[vv], self.Lk)
            for i in range(Lk):
                self.dkk[vv][i] = self.dkk[vv][i] / phase
        return self.dkk

    #normalize our results by the density of states of EPW calculation
    def NormVal(self, evals, dosef_k_dos, dosef_k_epc):
        dosef_single_spin_dos = sum(dosef_k_dos[:,0]) / 2 + sum(dosef_k_dos[:,1]) / 2
        dosef_single_spin_epc = sum(dosef_k_epc[:,0]) / 2 + sum(dosef_k_epc[:,1]) / 2
        evals = evals * dosef_single_spin_epc / dosef_single_spin_dos
        return evals
    
    def RunAll(self):
        
        #calculate the interaction matrix
        self.Vkk = self.CalVkk(self.Lk, self.Libnd, self.Vkk, self.Vkkdata1, self.dosef1, self.nmode)
        
        #setting how many eigenstates we calculate
        NE=36
        
        #calculate the eigenvector of interaction matrix Vkk
        evals, evecs = self.CalEigen(self.Vkk, NE)
        
        #Change the storage format of eigenvectors to store the four components
        #(singlet:\sigma_0; triplet:\sigma_1, \sigma_2, \sigma_3) of each k-point as a 2 * 2 matrix
        dkk = self.SaveToDelta(NE, evecs, self.Lk, self.Libnd)
        
        #because the define of density of states of dos calculation and epc calculation is different,
        #we normalize our results by the density of states of EPW calculation
        evals = self.NormVal(evals, self.dosef_dos, self.dosef1)
        
        #save our delta results
        self.SaveAll()

        print(evals)
#cal_delta(prefix, nk1,nk2, nk3, nmode)
new = cal_delta("pb", 20, 20, 20, 3)
new.RunAll()
