import ctypes

_br = ctypes.CDLL('libBAMreport.so')
_br.bamreports.argtypes = (ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_bool)

def bamreports(refName, alnName, rlnName, qulName, ncsName, covName, insName, gccName, stsName, windowsize):
    global _br
    _br.bamreports(ctypes.c_char_p(refName), 
                   ctypes.c_char_p(alnName),
                   ctypes.c_char_p(rlnName),
                   ctypes.c_char_p(qulName),
                   ctypes.c_char_p(ncsName),
                   ctypes.c_char_p(covName),
                   ctypes.c_char_p(insName),
                   ctypes.c_char_p(gccName),
                   ctypes.c_int(windowsize),
                   ctypes.c_char_p(stsName),
                   ctypes.c_char_p(None),
                   ctypes.c_char_p(None),
                   ctypes.c_bool(False))
                
    return
