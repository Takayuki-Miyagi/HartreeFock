import os
import sys

def get_nnfile(path, opr, pot, renorm, cut, hw, emax, e2max, txb, fom):
    f = path + '/' + opr + '-HO_NN-only_' + pot + '_' + renorm
    if(renorm != 'bare'):
        f += cut
    f += '_hw' + str(hw) + '_emax' + str(emax) + '_e2max' + str(e2max)
    f += '.' + txb + '.' + fom
    return f

def get_nnnfile(path, opr, renorm, cut, cd, ce, lam, e3max, hw, gen, txb):
    f = path + '/' + opr + '-' + renorm
    if(renorm != 'bare'):
        f += cut
    if(gen):
        f += '_cD' + cd + 'cE' + ce + '_lam' + lam
    f += '_e3max' + str(e3max) + '_hw' + str(hw)
    if(renorm != 'bare'):
        if(gen):
            f += '_NNN-full'
        else:
            f += '_NNN-ind'
    f += '.' + txb
    return f

def get_summary_file(Z, N, A, HFloop, pot, renorm, cut, gen, hw, emax, e2max, thbme):
    if(HFloop):
        f = 'Summary-HF-'
    else:
        f = 'Summary-HO-'
    f += 'Z' + str(Z) + 'N' + str(N) + 'A' + str(A) + '_' + pot + '_' + renorm
    if(renorm != 'bare'):
        f += cut
    if(thbme):
        if(gen):
            f += '_NN+3N'
        else:
            f += '_NN'
    else:
        f += '_NN-only'
    f += '_hw' + str(hw) + '_emax' + str(emax) + '_e2max' + str(e2max)
    f += '.dat'
    return f

def get_hamil_file(Z, N, A, HFloop, pot, renorm, cut, gen, hw, emax, e2max, txb, thbme):
    if(HFloop):
        f = 'Hamil-HF-'
    else:
        f = 'Hamil-HO-'
    f += 'Z' + str(Z) + 'N' + str(N) + 'A' + str(A) + '_' + pot + '_' + renorm
    if(renorm != 'bare'):
        f += cut
    if(thbme):
        if(gen):
            f += '_NN+3N'
        else:
            f += '_NN'
    else:
        f += '_NN-only'
    f += '_hw' + str(hw) + '_emax' + str(emax) + '_e2max' + str(e2max)
    f += '.snt.' + txb
    return f
