#!/usr/local/bin/python3
import sys
import os.path
import subprocess
from collections import OrderedDict
HOME = os.path.expanduser('~')

Params=OrderedDict()
exe = 'HartreeFock.exe'
inputf='hf_input.txt'
hwlist=[24]
elist=[4,6]
e3list = [6]
NuclList=['O16']

emax_nn = 8
e2max_nn = 16
emax_3n = 8
e2max_3n = 8
e3max_3n = 8

Params['optrs']='Rm2,Rp2,Rn2'
Params['optrs']=''

def SetHamil(hw, emax_nn, e2max_nn, emax_3n, e2max_3n, e3max_3n):
    path2 = HOME + '/MtxElmnt/2BME'
    path3 = HOME + '/MtxElmnt/3BME'
    file_nn_int = 'TwBME-HO_NN-only_N3LO_EM500_srg2.00_hw'
    file_nn_int += str(hw)+'_emax'+str(emax_nn)+'_e2max'+str(e2max_nn)+'.txt.me2j'
    file_3n_int = 'none'
    return path2, file_nn_int, path3, file_3n_int

def main():
    for Nucl in NuclList:
        for e in elist:
            for e3 in e3list:
                for hw in hwlist:
                    path2, file_nn, path3, file_3n = SetHamil(hw,emax_nn,e2max_nn,emax_3n,e2max_3n,e3max_3n)
                    file_nn_int = path2 + '/' + file_nn
                    file_3n_int = 'none'
                    if(file_3n != 'none'):
                        file_3n_int = path3 + '/' + file_3n

                    Params['Nucl'] = Nucl
                    Params['hw']=hw
                    Params['emax'] = e
                    Params['e3max'] = e3
                    Params['emax_nn'] = emax_nn
                    Params['e2max_nn'] = e2max_nn
                    Params['e2max_nn'] = e2max_nn
                    Params['emax_3n'] = emax_3n
                    Params['e2max_3n'] = e2max_3n
                    Params['e3max_3n'] = e3max_3n
                    Params['int_nn_file'] = file_nn_int Params['int_3n_file'] = file_3n_int
                    f = open(inputf,'w')
                    f.write('&input \n')
                    for key, value in Params.items():
                        if(isinstance(value,str)):
                            f.write(str(key)+'= "'+value+'" \n')
                        else:
                            f.write(str(key)+'= '+str(value)+' \n')
                    f.write("&end \n")
                    f.close()
                    cmd = exe + ' ' + inputf
                    subprocess.call(cmd,shell=True)

if(__name__ == '__main__'):
    main()


