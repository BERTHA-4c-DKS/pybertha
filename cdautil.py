import numpy

#def join_cmat(cmata, ndima, nocca, cmatb, ndimb, noccb, ndimab, noccab):
#    intial check
#    if (cmata.shape == (ndima,nocca)) :
#       print("check cmata: OK")
#    else: 
#       print("cmata : wrong dim")
#    if (cmatb.shape == (ndimb,noccb)) :
#       print("check cmatb: OK")
#    else: 
#       print("cmatb : wrong dim")
#    ntot = ndima + ndimb
#    totocc = nocca + noccb
#    if ( ntot == ndimab) : 
#      print("cmatjoin[0] OK")
#    else:
#      print("check cmatjoin[0]")
#    if ( totocc == noccab) :
#      print("cmatjoin[1] OK")
#    else:
#      print("check cmatjoin[1]")

def join_cmat(cmata, cmatb, ndimab):
    if (cmata.shape[0] % 2 == 0) :
       print("check cmata: OK")
    else: 
       print("cmata : wrong ndima")
    if (cmatb.shape[0] % 2 == 0) :
       print("check cmatb: OK")
    else: 
       print("cmatb : wrong ndimb")
    ndimtot =  cmata.shape[0] + cmatb.shape[0]
    if (ndimtot == ndimab):
       print("check cmat_join: OK")
    else:
       print("check cmat_join : wrong ndimtot")
    noccab =  cmata.shape[1] + cmatb.shape[1]
    cmat_join = numpy.zeros((ndimab,noccab),dtype=numpy.complex128)
    
    ndima = cmata.shape[0]
    nocca = cmata.shape[1]
    ndimb = cmatb.shape[0]
    noccb = cmatb.shape[1]
  
    cmat_join[0:ndima/2,0:nocca] = cmata[0:ndima/2,0:nocca]
    cmat_join[ndimab/2:ndimab/2+ndima/2,0:nocca] = cmata[ndima/2:,0:nocca]
    cmat_join[ndima/2:ndimab/2,nocca:] = cmatb[0:ndimb/2,0:noccb]
    cmat_join[ndimab/2+ndima/2: ,nocca:] = cmatb[ndimb/2:,0:noccb]
    
    return cmat_join
