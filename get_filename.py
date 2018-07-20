import os
import sys

def get_nnfile(path, opr, pot, renorm, cut, hw, emax, e2max, txb, fom):
    f = path + '/' + opr + '-HO_NN-only_' + pot + '_' + renorm
    if(renorm != 'bare'):
        f += cut
    f += '_hw' + hw + '_emax' + str(emax) + '_e2max' + str(e2max)
    f += '.' + txb + '.' + fom
    return f

def get_nnnfile(path, opr, renorm, cut, cd, ce, lam, e3max, hw, gen, txb):
    f = path + '/' + opr + '-' + renorm + cut
    if(gen):
        f += '_cD' + cd + 'cE' + ce + '_lam' + lam
    f += '_e3max' + str(e3max) + '_hw' + str(hw)
    if(gen):
        f += '_NNN-full'
    else:
        f += '_NNN-ind'
    f += '.' + txb
    return f

def get_summary_file(pot, renorm, cut, gen, hw, emax, e2max):
    f = 'Summary-' + pot + '_' + renorm
    if(renorm != 'bare'):
        f += cut
    if(gen):
        f += '_NN+3N'
    if(gen):
        f += '_NN'
    f += '_hw' + str(hw) + '_emax' + str(emax) + '_e2max' + str(e2max)
    f += '.dat'

def get_hamil_file(pot, renorm, cut, gen, hw, emax, e2max, txb):
    f = 'Hamil-' + pot + '_' + renorm
    if(renorm != 'bare'):
        f += cut
    if(gen):
        f += '_NN+3N'
    if(gen):
        f += '_NN'
    f += '_hw' + str(hw) + '_emax' + str(emax) + '_e2max' + str(e2max)
    f += '.snt.' + txb
