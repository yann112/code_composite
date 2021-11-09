
# coding: utf-8

from lib_calcul_matrices import *
import numpy as np


##########################    Matrices de rigiditees   ##############################################
####copier coller les matrices de rigidité et les donnees d entrees


PO =  15 #charge en Pa 
a = 2 #longueur de plaque en metre
x_1 = 1#point etudie en metre


#contrainte en Mpa
SIG1_rt = 1100
SIG1_rc = -250
SIG2_rt = 35
SIG2_rc = -120
SIG6_r = 50


#Matrices de rigidite en Gpa
Q_pi_0 =Matrix([
[117.23, 2.08, 0],
[  2.08, 5.94, 0],
[     0,    0, 2]])

Q_pi_2 =Matrix([
[5.94,   2.08, 0],
[2.08, 117.23, 0],
[   0,      0, 2]])

Q_pi_4 = Matrix([
[33.83, 29.83, 27.82],
[29.83, 33.83, 27.82],
[27.82, 27.82, 29.75]])

_Q_pi_4 = Matrix([
[ 33.83,  29.83, -27.82],
[ 29.83,  33.83, -27.82],
[-27.82, -27.82,  29.75]])

D = Matrix([
[746.13,  623.69, 445.16],
[623.69, 1645.73, 445.16],
[445.16,  445.16, 621.46]])

#forme de u3O_11 pour rotule/rotule chargement uniforme : -PO/2/D[0,0]*(x_1**2-a*x_1))
#forme de u3O_11 pour encastrement/encastrement chargement -q1*sin(pi*x_1/l) : PO/D[0,0]* (a/np.pi)**2*(np.sin(np.pi*x_1/a)-2/np.pi)
#forme des composantes de eps pour une plaque en appui simple sur contour et chargement q11*sin(pi*x_1/a)*sin(pi*x_2/a)
##k = (a/np.pi)**2 * PO/(D[0,0]+D[1,1]+2*D[0,1]+4*D[2,2])
##u3O_11 = k*-sin(np.pi*x_1/a)*sin(np.pi*x_2/a)
##u3O_22 = k*-sin(np.pi*x_1/a)*sin(np.pi*x_2/a)
##u3O_12 = k*cos(np.pi*x_1/a)*cos(np.pi*x_2/a)


u3O_11 = - PO/12/D[0,0]* (6*x_1**2-6*a*x_1+a**2)
u3O_22 = 0
u3O_12 = 0

#remplacer les info pour chaque couches exemple c0_inf0 = [hauteur,angle,matrice]
#remplacer les couches inutiles par  [-0.0,0,zeros(3,3)]

c1_inf0 = [-3.5,pi/4,Q_pi_4]
c1_sup0 = [-2.5,pi/4,Q_pi_4]

c2_inf0 = [-2.5,pi/2,Q_pi_2]
c2_sup0 = [-1.5,pi/2,Q_pi_2]

c3_inf0 = [-1.5,-pi/4,_Q_pi_4]
c3_sup0 = [-0.5,-pi/4,_Q_pi_4]

c4_inf0 = [-0.5,0,Q_pi_0]
c4_sup0 = [-0.5,0,Q_pi_0]

c5_inf0 = [0.5,-pi/4,_Q_pi_4]
c5_sup0 = [1.5,-pi/4,_Q_pi_4]

c6_inf0 = [1.5,pi/2,Q_pi_2]
c6_sup0 = [2.5,pi/2,Q_pi_2,]

c7_inf0 = [2.5,pi/4,Q_pi_4]
c7_sup0 = [3.5,pi/4,Q_pi_4]





















####Partie calcul a ne pas toucher

print()
eps = Matrix([-u3O_11, -u3O_22, -2*u3O_12])

print("La matrice des deformation est epsilon = X3 * ")
print(eps)
print()
print("les contraintes des couches dans leurs propres base sont : ")
print()
c1_inf = c1_inf0[0]*c1_inf0[2]*eps
c1_sup = c1_sup0[0]*c1_sup0[2]*eps
c2_inf = c2_inf0[0]*c2_inf0[2]*eps
c2_sup = c2_sup0[0]*c2_sup0[2]*eps
c3_inf = c3_inf0[0]*c3_inf0[2]*eps
c3_sup = c3_sup0[0]*c3_sup0[2]*eps
c4_inf = c4_inf0[0]*c4_inf0[2]*eps
c4_sup = c4_sup0[0]*c4_sup0[2]*eps
c5_inf = c5_inf0[0]*c5_inf0[2]*eps
c5_sup = c5_sup0[0]*c5_sup0[2]*eps
c6_inf = c6_inf0[0]*c6_inf0[2]*eps
c6_sup = c6_sup0[0]*c6_sup0[2]*eps
c7_inf = c7_inf0[0]*c7_inf0[2]*eps
c7_sup = c7_sup0[0]*c7_sup0[2]*eps


print("c1_inf")
print(c1_inf)
print("c1_sup")
print(c1_sup)
print()

print("c2_inf")
print(c2_inf)
print("c2_sup")
print(c2_sup)
print()

print("c3_inf")
print(c3_inf)
print("c3_sup")
print(c3_sup)
print()

print("c4_inf")
print(c4_inf)
print("c4_sup")
print(c4_sup)
print()

print("c5_inf")
print(c5_inf)
print("c5_sup")
print(c5_sup)
print()

print("c6_inf")
print(c6_inf)
print("c6_sup")
print(c6_sup)
print()

print("c7_inf")
print(c7_inf)
print("c7_sup")
print(c7_sup)
print()
print()
print()
print()
print()






print("on rammene les contraintes dans la base d orthotropie avec un changement de base")
print()

c1_inf_h = cdb_sigma(c1_inf,c1_inf0[1])
c1_sup_h = cdb_sigma(c1_sup,c1_sup0[1])
c2_inf_h = cdb_sigma(c2_inf,c2_inf0[1])
c2_sup_h = cdb_sigma(c2_sup,c2_sup0[1])
c3_inf_h = cdb_sigma(c3_inf,c3_inf0[1])
c3_sup_h = cdb_sigma(c3_sup,c3_sup0[1])
c4_inf_h = cdb_sigma(c4_inf,c4_inf0[1])
c4_sup_h = cdb_sigma(c4_sup,c3_sup0[1])
c5_inf_h = cdb_sigma(c5_inf,c5_inf0[1])
c5_sup_h = cdb_sigma(c5_sup,c5_sup0[1])
c6_inf_h = cdb_sigma(c6_inf,c6_inf0[1])
c6_sup_h = cdb_sigma(c6_sup,c6_sup0[1])
c7_inf_h = cdb_sigma(c7_inf,c7_inf0[1])
c7_sup_h = cdb_sigma(c7_sup,c7_sup0[1])


print("c1_inf")
print(c1_inf_h)
print("c1_sup")
print(c1_sup_h)
print()

print("c2_inf")
print(c2_inf_h)
print("c2_sup")
print(c2_sup_h)
print()

print("c3_inf")
print(c3_inf_h)
print("c3_sup")
print(c3_sup_h)
print()

print("c4_inf")
print(c4_inf_h)
print("c4_sup")
print(c4_sup_h)
print()

print("c5_inf")
print(c5_inf_h)
print("c5_sup")
print(c5_sup_h)
print()

print("c6_inf")
print(c6_inf_h)
print("c6_sup")
print(c6_sup_h)
print()

print("c7_inf")
print(c7_inf_h)
print("c7_sup")
print(c7_sup_h)
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
print("pour c1_inf")
cont_max(c1_inf_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c1_inf_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()
print("pour c1_sup")
cont_max(c1_sup_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c1_sup_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()
print()
print()
print("pour c2_inf")
cont_max(c2_inf_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c2_inf_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()
print("pour c2_sup")
cont_max(c2_sup_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c2_sup_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()
print()
print()
print("pour c3_inf")
cont_max(c3_inf_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c3_inf_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()
print("pour c3_sup")
cont_max(c3_sup_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c3_sup_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()
print()
print()
print("pour c4_inf")
cont_max(c4_inf_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c4_inf_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()
print("pour c4_sup")
cont_max(c4_sup_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c4_sup_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()
print()
print()
print("pour c5_inf")
cont_max(c5_inf_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c5_inf_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()
print("pour c5_sup")
cont_max(c5_sup_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c5_sup_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()
print()
print()
print("pour c6_inf")
cont_max(c6_inf_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c6_inf_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()
print("pour c6_sup")
cont_max(c6_sup_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c6_sup_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()
print()
print()
print("pour c7_inf")
cont_max(c7_inf_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c7_inf_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
print()
print("pour c7_sup")
cont_max(c7_sup_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
tsai_hill(c7_sup_h,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r )
      

pyplot.plot([c1_inf[0],c1_sup[0]], [c1_inf0[0],c1_sup0[0]])
pyplot.plot([c2_inf[0],c2_sup[0]], [c2_inf0[0],c2_sup0[0]])
pyplot.plot([c3_inf[0],c3_sup[0]], [c3_inf0[0],c3_sup0[0]])
pyplot.plot([c4_inf[0],c4_sup[0]], [c4_inf0[0],c4_sup0[0]])
pyplot.plot([c5_inf[0],c5_sup[0]], [c5_inf0[0],c5_sup0[0]])
pyplot.plot([c6_inf[0],c6_sup[0]], [c6_inf0[0],c6_sup0[0]])
pyplot.plot([c7_inf[0],c7_sup[0]], [c7_inf0[0],c7_sup0[0]])
pyplot.ylabel('hauteurs des couches')
pyplot.xlabel('contraintes')
pyplot.suptitle('contraintes sigma 1')
pyplot.show()
pyplot.close()

pyplot.plot([c1_inf[1],c1_sup[1]], [c1_inf0[0],c1_sup0[0]])
pyplot.plot([c2_inf[1],c2_sup[1]], [c2_inf0[0],c2_sup0[0]])
pyplot.plot([c3_inf[1],c3_sup[1]], [c3_inf0[0],c3_sup0[0]])
pyplot.plot([c4_inf[1],c4_sup[1]], [c4_inf0[0],c4_sup0[0]])
pyplot.plot([c5_inf[1],c5_sup[1]], [c5_inf0[0],c5_sup0[0]])
pyplot.plot([c6_inf[1],c6_sup[1]], [c6_inf0[0],c6_sup0[0]])
pyplot.plot([c7_inf[1],c7_sup[1]], [c7_inf0[0],c7_sup0[0]])
pyplot.ylabel('hauteurs des couches')
pyplot.xlabel('contraintes')
pyplot.suptitle('contraintes sigma 2')
pyplot.show()
pyplot.close()

pyplot.plot([c1_inf[2],c1_sup[2]], [c1_inf0[0],c1_sup0[0]])
pyplot.plot([c2_inf[2],c2_sup[2]], [c2_inf0[0],c2_sup0[0]])
pyplot.plot([c3_inf[2],c3_sup[2]], [c3_inf0[0],c3_sup0[0]])
pyplot.plot([c4_inf[2],c4_sup[2]], [c4_inf0[0],c4_sup0[0]])
pyplot.plot([c5_inf[2],c5_sup[2]], [c5_inf0[0],c5_sup0[0]])
pyplot.plot([c6_inf[2],c6_sup[2]], [c6_inf0[0],c6_sup0[0]])
pyplot.plot([c7_inf[2],c7_sup[2]], [c7_inf0[0],c7_sup0[0]])
pyplot.ylabel('hauteurs des couches')
pyplot.xlabel('contraintes')
pyplot.suptitle('contraintes sigma 6')
pyplot.show()
pyplot.close()
