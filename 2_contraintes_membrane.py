
# coding: utf-8

from lib_calcul_matrices import *
import numpy as np


##########################    Matrices de rigiditees   ##############################################
####copier coller les matrices de rigidité et les donnees d entrees

#contrainte en Mpa
SIG1_rt = 1200
SIG1_rc = -300
SIG2_rt = 50
SIG2_rc = -140
SIG6_r = 70


#Matrices de rigidite en Gpa
Q_pi_0 =Matrix([
[70.5, 2.01, 0],
[2.01, 8.06, 0],
[   0,    0, 5]])

Q_pi_2 = Matrix([
[8.06, 2.01, 0],
[2.01, 70.5, 0],
[   0,    0, 5]])

#Matrice de rigidite en membrane en N.m
A = Matrix([
[745323741.01,  30215827.34,        0],
[ 30215827.34, 433093525.18,        0],
[           0,            0, 75000000]])

#Matrice des efforts en N.m
N = Matrix([10**6,0,0]).reshape(1,3)

#remplacer les info pour chaque couches exemple c0_inf0 = [angle,matrice]
#remplacer les couches inutiles par  [-0.0,0,zeros(3,3)]

c1_inf0 = [0,zeros(3,3)]
c1_sup0 = [0,zeros(3,3)]

c2_inf0 = [0,zeros(3,3)]
c2_sup0 = [0,zeros(3,3)]

c3_inf0 = [0,Q_pi_0]
c3_sup0 = [0,Q_pi_0]

c4_inf0 = [pi/2,Q_pi_2]
c4_sup0 = [pi/2,Q_pi_2]

c5_inf0 = [0,Q_pi_0]
c5_sup0 = [0,Q_pi_0]

c6_inf0 = [0,zeros(3,3)]
c6_sup0 = [0,zeros(3,3)]

c7_inf0 = [0,zeros(3,3)]
c7_sup0 = [0,zeros(3,3)]





















####Partie calcul a ne pas toucher

print()
A_1 = A.inv()
print("L'inverse de la matrice A est : ")
pprint(A_1)
eps_0 = N*A_1
eps = eps_0.T
print()
print("On en déduit grace à la relation N = A*eps donc eps = N * A-1 la matrice des deformation epsilon * ")
pprint(eps)
print()
print("les contraintes des couches dans leurs propres base sont C * eps en Mpa: ")
print()
c1 = c1_inf0[1]*eps *10**3 #le 10**3 est pour avoir des Mpa
c2 = c2_inf0[1]*eps *10**3
c3 = c3_inf0[1]*eps *10**3
c4 = c4_inf0[1]*eps *10**3
c5 = c5_inf0[1]*eps *10**3
c6 = c6_inf0[1]*eps *10**3
c7 = c7_inf0[1]*eps *10**3



print("c1")
print(c1)

print()

print("c2")
print(c2)

print()

print("c3")
print(c3)

print()

print("c4")
print(c4)

print()

print("c5")
print(c5)

print()

print("c6")
print(c6)

print()

print("c7")
print(c7)

print()
print()
print()
print()
print()






print("on rammene les contraintes dans la base d orthotropie avec un changement de base")
print()

c1_h = cdb_sigma(c1,c1_inf0[0])
c2_h = cdb_sigma(c2,c2_inf0[0])
c3_h = cdb_sigma(c3,c3_inf0[0])
c4_h = cdb_sigma(c4,c4_inf0[0])
c5_h = cdb_sigma(c5,c5_inf0[0])
c6_h = cdb_sigma(c6,c6_inf0[0])
c7_h = cdb_sigma(c7,c7_inf0[0])



print("c1")
print(c1_h)

print()

print("c2")
print(c2_h)

print()

print("c3")
print(c3_h)

print()

print("c4")
print(c4_h)

print()

print("c5")
print(c5_h)

print()

print("c6")
print(c6_h)

print()

print("c7")
print(c7_h)

print()
print()
print()
print()
print()
print("Calcul des coefficients de chargements:")
print()
print()
print("pour les critères de ruptures suivants: sig1rt = %s sig1rc = %s sig2rt = %s sig2rc = %s sig6 = %s" % (SIG1_rt, SIG1_rc, SIG2_rt, SIG2_rc, SIG6_r))
print()
print("pour c1")
cont_max(c1_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c1_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )

print()
print()
print()
print("pour c2")
cont_max(c2_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c2_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()

print()
print()
print()
print("pour c3")
cont_max(c3_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c3_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )

print()
print()
print()
print("pour c4")
cont_max(c4_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c4_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()

print()
print()
print()
print("pour c5")
cont_max(c5_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c5_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()

print()
print()
print()
print("pour c6")
cont_max(c6_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c6_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()

print()
print()
print()
print("pour c7")
cont_max(c7_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c7_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()

      

