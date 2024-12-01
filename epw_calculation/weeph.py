#from builtins import input
import numpy as np
import os
import math
import numpy.linalg as la
ryd2ev=13.6
def is_empty_file_1(file_path: str):
    assert isinstance(file_path, str), f"file_path: {type(file_path)}"
    assert os.path.isfile(file_path), f"file_path: {file_path}"
    return os.path.getsize(file_path) == 0
def is_empty_file_4(file_path: str):
    assert isinstance(file_path, str), f"file_path: {type(file_path)}"
    assert os.path.isfile(file_path), f"file_path: {file_path}"

    with open(file_path, "r", encoding="utf-8") as f:
        first_char = f.read(1)
        if not first_char:
            return False
    with open(file_path, "r", encoding="utf-8") as f:
        first_char = f.read(1)
        if first_char:
            return True
#    return False
def get_unique_N(iterable):
    """Yields (in order) the first N unique elements of iterable. 
    Might yield less if data too short."""
    seen = set()
    for e in iterable:
        if e in seen:
            continue
        seen.add(e)
        yield e

#print(np.loadtxt('dos.dat'))
def fff():
    prefix="pb"
    here=os.getcwd()
    ephipool=here+'/'+prefix+'.ephmat'

    #merge all data from every .ephmat file
    os.chdir(ephipool)
    ii=0
    Nppp=64
    N=1+np.arange(Nppp)
    npool=0
    jj=1
    while jj<Nppp:
        print(jj)
        ipool=jj
        if is_empty_file_4('ephmat'+str(ipool)):
            if len(np.array(np.loadtxt('ephmat'+str(ipool))))>2:
                VV = np.array(np.loadtxt('ephmat'+str(ipool)))
                break
            else: jj=jj+1
        else:
            jj=jj+1
        
    while ipool<=Nppp:
        ii=ipool+1

        while ii<=Nppp:
            if is_empty_file_4('ephmat'+str(ii)):
                print(ii)
                print(is_empty_file_4('ephmat'+str(ii)))
                print(np.loadtxt('ephmat'+str(ii)))
                if len(np.array(np.loadtxt('ephmat'+str(ii))))>2:
                    break
                else:
                    ii=ii+1
            else:
                print(ii)
                ii=ii+1

        if ii>Nppp:
            break#i
        ipool=ii
        print(ii)
        print(is_empty_file_4('ephmat'+str(ipool)))

        ephmatfile=np.loadtxt('ephmat'+str(ipool))
        ff = np.array(np.loadtxt('ephmat'+str(ipool)))
        VV=np.concatenate((VV,ff),axis=0)


    #get data size of input file
    nbndsub=len(np.unique(VV[:,1]))
    nmode=6
    Lk=len(np.unique(VV[:,0]))

    #the index of klist is in fine k mesh but index of ktot is the same with the index of lambdaFS
    #The k-point order of lambda FS is obtained by sorting the k-points near the Fermi surface according to the index size in klist
    klist=np.array(np.unique(VV[:,0],True))
    ktot=np.zeros(int(max(klist[0])))
    for i in range(len(klist[0])):
        ktot[int(klist[0][i]-1)]=i
        #print(klist[0][i],i)

    #save wave function data
    wfc=np.zeros((Lk,nbndsub,nbndsub,nmode,2))

    #write wave function data
    for i in range(len(VV)):
        k=int(VV[i,0]-1)
        ik=int(ktot[k]-1)
        ibnd=int(VV[i,1]-1)
        jbnd=int(VV[i,2]-1)
        imode=int(VV[i,3]-1)
        #print(k,ik)
        wfc[ik,ibnd,jbnd,imode,0]=VV[i,4]
        wfc[ik,ibnd,jbnd,imode,1]=VV[i,5]
    os.chdir(here)
    f=open('wfcfff.dat','w')
    for i in range(len(VV)):
        k=int(VV[i,0]-1)
        ik=int(ktot[k])
        ibnd=int(VV[i,1]-1)
        jbnd=int(VV[i,2]-1)
        imode=int(VV[i,3]-1)
        #print(k,ik)
        wfc[ik,ibnd,jbnd,imode,0]=VV[i,4]
        wfc[ik,ibnd,jbnd,imode,1]=VV[i,5]
        f.write(str(ik)+'   '+str(ibnd+1)+'   '+str(jbnd+1)+'   '+str(imode+1)+'   '+str(wfc[ik,ibnd,jbnd,imode,0])+'   '+str(wfc[ik,ibnd,jbnd,imode,1])+'\n')
    os.chdir(here)
    f.close()



    f=open('wfc333.dat','w')
    for ik in range(Lk):
        for ibnd in range(nbndsub):
            for imode in range(nmode):
                for jbnd in range(nbndsub):
                    f.write(str(ik+1)+'   '+str(ibnd+1)+'   '+str(jbnd+1)+'   '+str(imode+1)+'   '+str(wfc[ik,ibnd,jbnd,imode,0])+'   '+str(wfc[ik,ibnd,jbnd,imode,1])+'\n')
    #print(klist[0])
    f.close()
    
fff()
