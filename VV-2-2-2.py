#!~/.conda/envs/nlp/bin/python3.7
import numpy as np
import matplotlib.pyplot as plt
import math as mmm
import numpy.linalg as la
from scipy.sparse.linalg.eigen.arpack import eigs as largest_eigsh
import csv

# import pandas as pd
# It then extracts the band data, and plots the bands, the fermi energy in red, and the high symmetry points
plt.rc("font", family='Times New Roman')

import os
prefix='Pb'
nmode = 3
nk1=20
nk2=20
nk3=20


Vkkfile = prefix+'.lambda_kkq'
Vkkdata = np.loadtxt(Vkkfile)
LVkk = len(Vkkdata)
kkkk = np.unique(Vkkdata[:, 0])

aaaa = np.unique(Vkkdata[:, 4])



datafile = prefix+'.lambda_FS'
lambdaFS = np.loadtxt(datafile)
kindex = np.unique(lambdaFS[:, 3], return_index=True)
Lk = len(kindex[0])
Libnd = len(aaaa)
Lkk = np.zeros((Lk))
kpoints = lambdaFS[kindex[1], 0:3]  # This is all the unique x-points
print(Lk)
save_path = 'delta'

Vkkdata1 = np.zeros((Lk, Lk, Libnd, Libnd, nmode), dtype=np.complex64)
Vkkdata2 = np.zeros((Lk, Lk, Libnd, Libnd, Libnd, Libnd), dtype=np.complex64)

dosef0 = np.zeros((Lk, Libnd))
dosef = np.zeros((Lk, Libnd))
dosef1 = np.zeros((Lk, Libnd))

Vkk = np.zeros((Lk * Libnd * Libnd, Lk * Libnd * Libnd), dtype=np.complex64)
Vkkss = np.zeros((Lk * Libnd * Libnd, Lk * Libnd * Libnd), dtype=np.complex64)

dkk = np.zeros((Lk, Libnd, Libnd), dtype=np.complex64)
orderk1 = np.zeros(LVkk)
orderk2 = np.zeros(LVkk)

#nk1 = 5
#nk2 = 5
#nk3 = 5

sum1 = 0


def totuple(a):
    try:
        return tuple(totuple(i) for i in a)
    except TypeError:
        return a


dos = 0
dodoseff0 = 0
dodoseff1 = 0

def writedos(LV, npool, iii):
    # sum1=0
    istart = iii * int(LV / npool)
    if iii < npool - 1:
        iend = istart + int(LV / npool)
    if iii == npool - 1:
        iend = LV - istart
    ii = istart
    while ii < iend:
        sk1 = int(Vkkdata[ii, 4]) - 1
        skq1 = int(Vkkdata[ii, 5]) - 1
        ik = int(Vkkdata[ii, 0]) - 1
        ikq = int(Vkkdata[ii, 2]) - 1

        dosef0[ik, sk1] = 1 / nk2 / nk2 / nk3 * Vkkdata[ii, 10]
        ii = ii + 1
    ii = istart
    while ii < iend:
        sk1 = int(Vkkdata[ii, 4]) - 1
        skq1 = int(Vkkdata[ii, 5]) - 1
        ik = int(Vkkdata[ii, 0]) - 1
        ikq = int(Vkkdata[ii, 2]) - 1
        iskq10 = skq1 - 2 * int(skq1 / 2)
        iskq11 = 2 * int(skq1 / 2)
        # dosef[ik,sk1]=1/nk2/nk2/nk3*Vkkdata[ii,13]
        dosikq = dosef0[ikq, iskq11 + iskq10] / 2 + dosef0[ikq, iskq11 + 1 - iskq10] / 2
        dosef1[ikq, skq1] = dosikq
        dosef[ik,sk1]=1/nk2/nk2/nk3*Vkkdata[ii,12]
        # dosef[ikq,iskq11+1-iskq10]=1/nk2/nk2/nk3*Vkkdata[ii,14]
        #dosef = Vkkdata[ii, 10]
        ii = ii + 1

def writeVkkdata(LV, npool, iii):
    sum1 = 0
    dos = 0
    istart = iii * int(LV / npool)
    if iii < npool - 1:
        iend = istart + int(LV / npool)
    if iii == npool - 1:
        iend = LV - istart
    ii = istart
    f=open('gkk.txt','w')
    while ii < iend:
        sk1 = int(Vkkdata[ii, 4]) - 1
        skq1 = int(Vkkdata[ii, 5]) - 1
        # print(Vkkdata[ii,8])
        imode = int(Vkkdata[ii, 8]) - 1
        ik = int(Vkkdata[ii, 0]) - 1
        ikq = int(Vkkdata[ii, 2]) - 1
        g1 = Vkkdata[ii, 6] + 1j * Vkkdata[ii, 7]
        dos = dos + 1 / nk2 / nk2 / nk3 * Vkkdata[ii, 10]
        if Vkkdata[ii, 9] > 0.0000001:
            # print(ik,ikq,sk1,skq1,imode)
            Vkkdata1[ik, ikq, sk1, skq1, imode] = g1 / np.sqrt(Vkkdata[ii, 9])
            sum1 = sum1 + abs(g1 * g1 / (nk1 * nk2 * nk3) / (nk1 * nk2 * nk3) * Vkkdata[ii, 10] * Vkkdata[ii, 11]) / \
                   Vkkdata[ii, 9]
            # print(sum1)
        kx=lambdaFS[int(kindex[1][ik]), 0]
        ky=lambdaFS[int(kindex[1][ik]),1]
        kz=lambdaFS[int(kindex[1][ik]),2]
        kqx=lambdaFS[int(kindex[1][ikq]),0]
        kqy=lambdaFS[int(kindex[1][ikq]),1]
        kqz=lambdaFS[int(kindex[1][ikq]),2]
        if kx>0.5:
            kx=kx-0.5
        if kqx>0.5:
            kqx=kqx-0.5
        if ky>0.5:
            ky=ky-0.5
        if kz>0.5:
            kz=kz-0.5
        if kqy>0.5:
            kqy=kqy-0.5
        if kqz>0.5:
            kqz=kqz-0.5
        llk=np.sqrt((kx-kqx)*(kx-kqx)+(ky-kqy)*(ky-kqy)+(kz-kqz)*(kz-kqz))
        if llk<3*1/nk1 and imode==2:
            f.write(('{:10.6f}'*6).format(kx,ky,kz,kqx,kqy,kqz))
            f.write(('{:10d}'*4).format(ik,ikq,sk1,skq1))
            f.write(('{:10.6f}'*2).format(Vkkdata[ii,6],Vkkdata[ii,7]))
            #print('\n')
            f.write('\n')
        ii = ii + 1
    f.close()
    print(sum1, dos, sum1 / np.sum(dosef) * 2)


iend = 0




def getphase(dkk):
    for i in range(Lk):
        if abs(dkk[i][0, 0]) > 0:
            a = dkk[i][0, 0] / abs(dkk[i][0, 0])
            #print(abs(dkk[i][0, 0]))
            break
    return a


def writedelta(dkk, aaaaa):
    eev = np.zeros((Lk * Libnd * Libnd), dtype=np.complex64)

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
                    #eev[i * Libnd * Libnd + ii*4+ibnd*2 + jbnd] = newdeltakk[ibnd, jbnd]
                    eev[i * Libnd * Libnd + ii*4+ibnd*2 + jbnd] = dkk[i][ibnd, jbnd]
            #print(ii)
            #print(newdeltakk)
    np.savetxt(save_path + '/' + 'dkk' + str(aaaaa) + '.txt', np.column_stack([eev.real, eev.imag]))
def writegkk(dkk, aaaaa):
    eev = np.zeros((Lk * Libnd * Libnd), dtype=np.complex64)

    for i in range(Lk):
        for ii in range(int(Libnd/2)):
            ik=kindex[1][i]
            kk1 = int(lambdaFS[ik, 3]) - 1
            kkk1 = int(lambdaFS[ik, 4]) - 1
            
            deltakk = np.array([[dkk[kk1][2*ii, 2*ii], dkk[kk1][2*ii, 2*ii+1]], [dkk[kk1][2*ii+1, 2*ii], dkk[kk1][2*ii+1, 2*ii+1]]])
            deltakk=deltakk+np.transpose(np.conjugate(deltakk))
            #print(deltakk)
            [eigkk, diakk] = np.linalg.eig(deltakk)

            if np.linalg.det(diakk) == -1:
                signkk = -1
            else:
                signkk = 1
            newdeltakk = signkk * np.dot(np.dot(np.conjugate(np.transpose(diakk)), deltakk), diakk)




            deltakkk = np.array([[dkk[kkk1][2*ii, 2*ii], dkk[kkk1][2*ii, 2*ii+1]], [dkk[kkk1][2*ii+1, 2*ii], dkk[kkk1][2*ii+1, 2*ii+1]]])
            deltakkk=deltakkk+np.transpose(np.conjugate(deltakkk))
            #print(deltakkk)
            [eigkkk, diakkk] = np.linalg.eig(deltakkk)

            if np.linalg.det(diakkk) == -1:
                signkkk = -1
            else:
                signkkk = 1
            newdeltakkk = signkkk * np.dot(np.dot(np.conjugate(np.transpose(diakkk)), deltakkk), diakkk)


            a = 0
            for imode in range(nmode):
                if abs(Vkkdata1[kk1, kkk1, 0, 0, imode]) + abs(Vkkdata1[kk1, kkk1, 0, 1, imode]) >= a:
                    a = abs(Vkkdata1[kk1, kkk1, 0, 0, imode]) + abs(Vkkdata1[kk1, kkk1, 0, 1, imode])
                    amode = imode
            matrix = [[Vkkdata1[kk1, kkk1, 2*ii, 2*ii, amode], Vkkdata1[kk1, kkk1, 2*ii, 2*ii+1, amode]],[Vkkdata1[kk1, kkk1, 2*ii+1, 2*ii, amode], Vkkdata1[kk1, kkk1, 2*ii+1, 2*ii+1, amode]]]
            diagkk = -signkk * np.dot(np.dot(np.conjugate(np.transpose(diakk)), matrix), diakkk)



            for ibnd in range(2):
                for jbnd in range(2):
                    eev[i * Libnd * Libnd + ii*4+ibnd*2 + jbnd] = diagkk[ibnd, jbnd]
            print(kk1,kkk1)
            print(matrix)
            print(diagkk)
    np.savetxt(save_path + '/' + 'gkk' + str(aaaaa) + '.txt', np.column_stack([eev.real, eev.imag]))


def bndplot(labels):
    # sum1=0
    # sum2=0

    ikk1 = 0
    for i in range(len(lambdaFS)):
        kk1 = int(lambdaFS[i, 3]) - 1
        kkk1 = int(lambdaFS[i, 4]) - 1
        if kk1 != kkk1:
            ikk1 = ikk1 + 1
        if ikk1 > 20:
            break

    #print(kk1, kkk1)
    print(Lk, Libnd)

    Lk1 = 0
    fthick = 0.4

    ii = 0

    ii = 0
    V00 = 0
    V01 = 0
    V10 = 0
    V11 = 0
    writedos(LVkk, 1, 0)
    writeVkkdata(LVkk, 1, 0)

    d1 = np.sum(dosef0[:, 0]) / 2
    d2 = np.sum(dosef0[:, 1]) / 2

    # for ik in range(Lk):
    #    for ikq in range(Lk):
    #        if abs(Vkkdata1[ik,ikq,0,0,0,0])<abs(Vkkdata1[ik,ikq,0,0,1,1]):
    #            for sk1 in range(2):
    #                for sk2 in range(2):
    #                    for skq1 in range(2):
    #                        for skq2 in range(2):
    #                            Vkkdata1[ik,ikq,sk1,sk2,skq1,skq2]=Vkkdata2[ik,ikq,sk1,sk2,1-skq1,1-skq2]
    dodosef0 = 0
    dodosef1 = 0
    for ik in range(Lk):
        dodosef0 = dodosef0 + dosef[ik]
        dodosef1 = dodosef1 + dosef1[ik]
        # if Vkkdata1[ik,ikq,sk1,sk2,skq1,skq2]!=0:
        #    Lk1=Lk1+1
        for ikq in range(Lk):
            for sk1 in range(Libnd):
                for sk2 in range(Libnd):
                    for skq1 in range(Libnd):
                        for skq2 in range(Libnd):
                            # if kpoints[ik,0]!=0 and kpoints[ik,1]!=0 and kpoints[ik,2]!=0 and kpoints[ikq,0]!=0 and kpoints[ikq,1]!=0 and kpoints[ikq,2]!=0:
                            # if sk1!=skq1 and sk2!=skq2:
                            #    Vkkdata1[ik,ikq,sk1,sk2,skq1,skq2]=Vkkdata1[ik,ikq,sk1,sk2,skq1,skq2]
                            # if sk1==skq1 and sk2==skq2:
                            ikk = int(lambdaFS[int(kindex[1][ik]), 4] - 1)
                            ikkq = int(lambdaFS[int(kindex[1][ikq]), 4] - 1)
                            for imode in range(nmode):
                                #if Vkkdata1[ik,ikq,sk1,skq1,imode]!=0 and Vkkdata1[ik,ikq,2*int(sk1/2)+(1-sk1%2),2*int(skq1/2)+(1-skq1%2),imode]==0:
                                #    Vkkdata1[ik,ikq,sk1,skq1,imode]=0
                                Vkk[ik*Libnd*Libnd+sk1*Libnd+sk2,ikq*Libnd*Libnd+(skq1)*Libnd+skq2]=Vkk[ik*Libnd*Libnd+sk1*Libnd+sk2,ikq*Libnd*Libnd+(skq1)*Libnd+skq2]-Vkkdata1[ik,ikq,sk1,skq1,imode]*np.conjugate(Vkkdata1[ik,ikq,sk2,skq2,imode])*dosef1[ikq,skq1]*2
                            phase=0


    # print(Vkk)
    ii = 0

    for i in range(len(Vkk)):
        for j in range(len(Vkk)):
            Vkkss[i, j] = Vkk[i, j]
    ii = 0
    ii = 0

    vvv = 0

    while ii < len(Vkk[:, 0]):
        jj = 0
        while jj < len(Vkk[:, 0]):
            # print(Vkk[ii,jj]-np.conjugate(Vkk[jj,ii]))
            jj = jj + 1
        ii = ii + 1
    # print(Vkk-np.transpose(np.conjugate(Vkk)))
    # print(Vkk)
    # Vkk=[[0,1-1j],[1+1j,0]]
    NE = 36
    eaa, evv = largest_eigsh(Vkkss, NE, which='SR')
    evals, evecs = largest_eigsh(Vkkss, NE, which='SR')
    #evals, evecs = la.eig(Vkkss)
    #print(evals)
    # print(evecs)
    # 把对应的gap-function存入文件

    # 建立存function的文件夹
    if not os.path.exists(save_path):
        os.mkdir(save_path)

    for vv in range(NE):
        for i in range(Lk):
            for iibnd in range(Libnd):
                for jjbnd in range(Libnd):
                    dkk[i][iibnd, jjbnd] = evecs[i * Libnd * Libnd + iibnd * Libnd + jjbnd, vv]
        phase = getphase(dkk)
        for i in range(Lk):
            dkk[i] = dkk[i] / phase
            #print(dkk[i])
        writedelta(dkk, vv)
        writegkk(dkk,vv)


    eevvv = np.loadtxt(save_path + '/' + 'dkk' + str(0) + '.txt').view(complex).reshape(-1)
    for i in range(Lk):
        for iii in range(int(Libnd/2)):
            for iibnd in range(2):
                for jjbnd in range(2):
                    dkk[i][2*iii+iibnd, 2*iii+jjbnd] = eevvv[i * Libnd * Libnd +iii*4+ iibnd * 2 + jjbnd]
        #print(dkk[i])

    aa0 = 0
    aa00 = 0
    aa1 = 0
    kline = np.linspace(0, len(dkk[:, 0, 0]), len(dkk[:, 0, 0]))
    # print(len(kline))
    # print(len(dkk[:,0,0]))
    plt.scatter(kline, abs(dkk[:, 0, 0]))
    # plt.title('Energy Spectrum of Chain with {} Sites'.format(Nsites))
    plt.show()
    plt.savefig('s-wave.png')
    #print(dkk[0])
    #print(dkk[1])
    vvv = 0

    for vv in range(NE):
        phase = 1
        for i in range(Lk):
            for iibnd in range(Libnd):
                for jjbnd in range(Libnd):
                    dkk[i][iibnd, jjbnd] = evecs[i * Libnd * Libnd + iibnd * Libnd + jjbnd, vv]
        phase = getphase(dkk)

        for i in range(Lk):
            dkk[i] = dkk[i] / phase
        # print(dkk[i])
        # print(vv,evals[vv])
        # print(dkk[kk1])
        # print(dkk[kkk1])
        for i in range(len(lambdaFS)):
            kk1 = int(lambdaFS[i, 3]) - 1
            kkk1 = int(lambdaFS[i, 4]) - 1
            if kk1 != kkk1 and abs(dkk[kk1][0, 0]) > 0.000003:
                break
        deltakk = np.array([[dkk[kk1][0, 0], dkk[kk1][0, 1]], [np.conjugate(dkk[kk1][0, 1]), dkk[kk1][1, 1]]])
        deltakkk = np.array([[dkk[kkk1][0, 0], dkk[kkk1][0, 1]], [np.conjugate(dkk[kkk1][0, 1]), dkk[kkk1][1, 1]]])
        #print(deltakk)
        #print(deltakkk)
        [eigkk, diakk] = np.linalg.eig(deltakk)
        [eigkkk, diakkk] = np.linalg.eig(deltakkk)

        if np.linalg.det(diakk) == -1:
            signkk = -1
        else:
            signkk = 1
        newdeltakk = -signkk * np.dot(np.dot(np.conjugate(np.transpose(diakk)), deltakk), diakk)
        if np.linalg.det(diakkk) == -1:
            signkkk = -1
        else:
            signkkk = 1
        newdeltakkk = -signkkk * np.dot(np.dot(np.conjugate(np.transpose(diakkk)), deltakkk), diakkk)
        print(vv, evals[vv])
        print(newdeltakk)
        print(newdeltakkk)

        a = 0
        for imode in range(nmode):
            if abs(Vkkdata1[kk1, kkk1, 0, 0, imode]) + abs(Vkkdata1[kk1, kkk1, 0, 1, imode]) >= a:
                a = abs(Vkkdata1[kk1, kkk1, 0, 0, imode]) + abs(Vkkdata1[kk1, kkk1, 0, 1, imode])
                amode = imode
        matrix = [[Vkkdata1[kk1, kkk1, 0, 0, amode], Vkkdata1[kk1, kkk1, 0, 1, amode]],
                  [Vkkdata1[kk1, kkk1, 1, 0, amode], Vkkdata1[kk1, kkk1, 1, 1, amode]]]

        # [eigkk,diakk]=np.linalg.eig(matrix)
        newdeltakk = -signkk * np.dot(np.dot(np.conjugate(np.transpose(diakk)), matrix), diakkk)
        print('gkk')
        print(matrix)
        print(abs(newdeltakk))
        kx = lambdaFS[kindex[1][i], 0]
        ky = lambdaFS[kindex[1][i], 1]
        kz = lambdaFS[kindex[1][i], 2]

    dodoseff0 = dodosef0[0] / 2 + dodosef0[1] / 2
    dodoseff1 = dodosef1[0] / 2 + dodosef1[1] / 2
    evals = evals * dodoseff1 / dodoseff0

    print(evals)
    #print(dkk[kk1])
    #print(dkk[kkk1])

    if np.linalg.det(diakk) == -1:
        signkk = -1
    else:
        signkk = 1

    # eq=np.unique(Vkkdata[:,15])
    # eq.sort()
    Lk1 = sum(Lkk)
    #print(Lk1)
    # dodoseff=dodosef[0]/2+dodosef[1]/2

    #print(V00 / d1, V11 / d2)
    #print(V01 / d1, V10 / d2)


bndplot("linewidth.png")
writeVkkdata(LVkk, 1, 0)
#print(sum1 / dodoseff0, sum1 / dodoseff0, sum1 / dodoseff1, sum2)
#print(dodoseff0, dodoseff0, dodoseff1)