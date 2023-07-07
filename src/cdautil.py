import numpy

def join_cmat(cmata, cmatb, ndimab):

    if (cmata.shape[0] % 2 != 0):
        return None

    if (cmatb.shape[0] % 2 != 0):
        return None

    ndimtot =  cmata.shape[0] + cmatb.shape[0]

    if (ndimtot != ndimab):
        return None

    noccab =  cmata.shape[1] + cmatb.shape[1]
    #print(("noccab : %i\n" % noccab))
    cmat_join = numpy.zeros((ndimab,noccab),dtype=numpy.complex128)
    
    ndima = cmata.shape[0]
    nocca = cmata.shape[1]
    ndimb = cmatb.shape[0]
    noccb = cmatb.shape[1]
  
    cmat_join[0:int(ndima/2),0:nocca] = cmata[0:int(ndima/2),0:nocca]
    cmat_join[int(ndimab/2):int(ndimab/2+ndima/2),0:nocca] = cmata[int(ndima/2):,0:nocca]
    cmat_join[int(ndima/2):int(ndimab/2),nocca:] = cmatb[0:int(ndimb/2),0:noccb]
    cmat_join[int(ndimab/2+ndima/2): ,nocca:] = cmatb[int(ndimb/2):,0:noccb]
    
    return cmat_join
