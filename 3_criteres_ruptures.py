# coding: utf-8

from lib_calcul_matrices import *


###Contraintes a rupture rt = traction rc = compression r = cisaillement en MPa
SIG1_rt = 1000
SIG1_rc = -800
SIG2_rt = 30
SIG2_rc = -100
SIG6_r = 70

###matrice de rigidite ramene dans la base d orthotropie en Mpa *10**3 si en Gpa
c0 =Matrix([[-33.2], [-9.38], [15.9]])
c1 = cdb_sigma(c0,-pi/6)
pprint(c1)
###critere de rupture

cont_max(c1,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c1,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
