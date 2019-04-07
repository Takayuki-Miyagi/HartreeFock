#!/usr/local/bin/python3
import sys
import os.path
import subprocess
from collections import OrderedDict
HOME = os.path.expanduser('~')

Params=OrderedDict()
exe = 'HartreeFock.exe'
inputf='hf_input.txt'
zlist=[1]
elist=[6]
Params['is_Atomic']=True
Params['is_MBPTEnergy']=False

def SetHamil(hw, emax_nn, e2max_nn, emax_3n, e2max_3n, e3max_3n):
    path2 = HOME + '/MtxElmnt/2BME'
    file_nn_int = 'AtomicHamil_LO_zeta'
    file_nn_int += str(hw)+'_emax'+str(emax_nn)+'_e2max'+str(e2max_nn)+'.snt'
    file_3n_int = 'none'
    return path2, file_nn_int, path3, file_3n_int

def main():
    for e in elist:
        for hw in hwlist:
            Params['Nucl'] = Nucl
            Params['hw']=hw
            Params['emax'] = e
            Params['e3max'] = e3
            Params['emax_nn'] = e
            Params['e2max_nn'] = 2*e
            path2, file_nn, path3, file_3n = SetHamil(hw,e,2*e,0,0,0)
            file_nn_int = path2 + '/' + file_nn
            file_3n_int = 'none'
            if(file_3n != 'none'):
                file_3n_int = path3 + '/' + file_3n
            Params['int_nn_file'] = file_nn_int
            Params['int_3n_file'] = file_3n_int
            f = open(inputf,'w')
            f.write('&input \n')
            for key, value in Params.items():
                if(isinstance(value,str)):
                    f.write(str(key)+'= "'+value+'" \n')
                else:
                    f.write(str(key)+'= '+str(value)+' \n')
            f.write("&end \n")
            f.close()
            cmd = exe + ' ' + inputf + " conf_atom_input.txt"
            subprocess.call(cmd,shell=True)
            f = "summary_"+file_nn
            cmd = "mv summary.out "+f
            subprocess.call(cmd,shell=True)


if(__name__ == '__main__'):
    main()


