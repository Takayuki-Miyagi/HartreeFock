class Orbit:
    def __init__(self):
        self.nlj_list = [{"n": 0, "l": 0, "j": 1},
                        {"n": 0, "l": 1, "j":  3},
                        {"n": 0, "l": 1, "j":  1},
                        {"n": 0, "l": 2, "j":  5},
                        {"n": 1, "l": 0, "j":  1},
                        {"n": 0, "l": 2, "j":  3},
                        {"n": 0, "l": 3, "j":  7},
                        {"n": 1, "l": 1, "j":  3},
                        {"n": 1, "l": 1, "j":  1},
                        {"n": 0, "l": 3, "j":  5},
                        {"n": 0, "l": 4, "j":  9},
                        {"n": 1, "l": 2, "j":  5},
                        {"n": 2, "l": 0, "j":  1},
                        {"n": 1, "l": 2, "j":  3},
                        {"n": 0, "l": 4, "j":  7},
                        {"n": 0, "l": 5, "j": 11},
                        {"n": 1, "l": 3, "j":  7},
                        {"n": 2, "l": 1, "j":  3},
                        {"n": 2, "l": 1, "j":  1},
                        {"n": 1, "l": 3, "j":  5},
                        {"n": 0, "l": 5, "j":  9},
                        {"n": 0, "l": 6, "j": 13},
                        {"n": 1, "l": 4, "j":  9},
                        {"n": 2, "l": 2, "j":  5},
                        {"n": 3, "l": 0, "j":  1},
                        {"n": 2, "l": 2, "j":  3},
                        {"n": 1, "l": 4, "j":  7},
                        {"n": 0, "l": 6, "j": 11},
                        ]

    def get_core(self,N,Z):
        core = []
        num = 0
        for nlj in self.nlj_list:
            n = nlj['n']
            l = nlj['l']
            j = nlj['j']
            num += j + 1
            if(num > N): continue
            nljtz = (n,l,j,1)
            core.append(nljtz)
        num = 0
        for nlj in self.nlj_list:
            n = nlj['n']
            l = nlj['l']
            j = nlj['j']
            num += j + 1
            if(num > Z): continue
            nljtz = (n,l,j,-1)
            core.append(nljtz)
        return core
    def get_valence(self,Nshell,Zshell):
        valence = []
        for nlj in self.nlj_list:
            n = nlj['n']
            l = nlj['l']
            j = nlj['j']
            if(2 * n + l != Nshell): continue
            nljtz = (n,l,j,1)
            valence.append(nljtz)
        for nlj in self.nlj_list:
            n = nlj['n']
            l = nlj['l']
            j = nlj['j']
            if(2 * n + l != Zshell): continue
            nljtz = (n,l,j,-1)
            valence.append(nljtz)
        return valence
    def get_no(self,N,Z,NO=None):
        no = {}
        num = 0
        for nlj in self.nlj_list:
            n = nlj['n']
            l = nlj['l']
            j = nlj['j']
            num += j + 1
            frac = 1
            if(num > N):
                if(NO != 'ENO'): break
                v = N - num + j + 1
                if(v > 0):
                    frac = float(v) / float(j+1)
                else:
                    frac = 0
                    break
            nljtz = (n,l,j,1)
            no[nljtz] = frac
        num = 0
        for nlj in self.nlj_list:
            n = nlj['n']
            l = nlj['l']
            j = nlj['j']
            num += j + 1
            frac = 1
            if(num > Z):
                if(NO != 'ENO'): break
                v = Z - num + j + 1
                if(v > 0):
                    frac = float(v) / float(j+1)
                else:
                    frac = 0
                    break
            nljtz = (n,l,j,-1)
            no[nljtz] = frac
        return no
