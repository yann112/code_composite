
# coding: utf-8
'''
test git
changement branche maison
'''

from lib_calcul_matrices import *

###########################Partie donnee d'entree a modifier selon le stratifie

####grandeur generalisée composite

E1 = 116.5#GPa
E2 = 5.9#GPa
v12 = 0.35#
G12 = 2#GPa

####grandeur generalisée ame en materiau homogene E1=E2 G = E/(2*(1+v))

##E = 5#GPa
##E = 5#GPa
##v = 0.25#
##G = E/(2*(1+v))#GPa  

####epaisseurs des couches
####inscrire 0 dans les couches non utilisee la centrale et les exterieures
e= 1

e1 = 1 
e2 = 1 
e3 = 1
e4 = 1
e5 = 1
e6 = 1 
e7 = 1 

####matrices de rigidité base d'orthotropie

c_base_ortho = cbarre(E1,E2,v12,G12)


####matrices de rigidité hors base d'orthotropie
####inscrire Matrix([0,0,0,0,0,0,0,0,0]).reshape(3,3) dans les couches non utilisee la centrale et les exterieures

c1 = ccouche("pi/4",pi/4,c_base_ortho)
c2 = ccouche("pi/2",pi/2,c_base_ortho)
c3 = ccouche("-pi/4",-pi/4,c_base_ortho)
c4 = ccouche("0",0,c_base_ortho)
c5 = c3
c6 = c2
c7 = c1







###########################Partie calcul   a ne pas toucher



####Matrices A

print()
print("Matrice de rigidité globale de membrane est A = en N.m-1")
print()
A = 10**6 * (e1*c1 + e2*c2 + e3*c3 + e4*c4 + e5*c5 + e6*c6 + e7*c7)
pprint(round_expr(A,2))
print()
print("soit 10^6 * ")
pprint(round_expr(A,2)/10**6)

####Matrices B

####hauteurs en mm  du centre de la couche /plan moyen
z4 = 0

z5 = z4 + 0.5*(e5+e4)
z6 = z5 + 0.5*(e6+e5)
z7 = z6 + 0.5*(e7+e6)

z3 = z4 - 0.5*(e4+e3)
z2 = z3 - 0.5*(e3+e2)
z1 = z2 - 0.5*(e2+e1)
print()

print()
print("Matrice de couplage membrane/flexion est B = en 10*3 Newtons")
B = c1*e1*z1 + c2*e2*z2 + c3*e3*z3 + c4*e4*z4 + c5*e5*z5 + c6*e6*z6 + c7*e7*z7
pprint(round_expr(B,2))



####Matrices D

print()
print("la suite des hauteurs des plans moyens est (attention aux couches virtuelles) :")
print(z1,z2,z3,z4,z5,z6,z7)

print()
print("Matrice de rigidité en flexion est D = en N.m")
D = 1/3 * ((c1*((z1+0.5*e1)**3 - (z1-0.5*e1)**3)) +
           (c2*((z2+0.5*e2)**3 - (z2-0.5*e2)**3)) +
           (c3*((z3+0.5*e3)**3 - (z3-0.5*e3)**3)) +
           (c4*((z4+0.5*e4)**3 - (z4-0.5*e4)**3)) +
           (c5*((z5+0.5*e5)**3 - (z5-0.5*e5)**3)) +
           (c6*((z6+0.5*e6)**3 - (z6-0.5*e6)**3)) +
           (c7*((z7+0.5*e7)**3 - (z7-0.5*e7)**3)))

pprint(round_expr(D,2))

####plot stratifie

list_h = (z1+0.5*e1 , z1-0.5*e1 ,
          z2+0.5*e2 , z2-0.5*e2 ,
          z3+0.5*e3 , z3-0.5*e3 ,
          z4+0.5*e4 , z4-0.5*e4 ,
          z5+0.5*e5 , z5-0.5*e5 ,
          z6+0.5*e6 , z6-0.5*e6 ,
          z7+0.5*e7 , z7-0.5*e7 ,)
for i in list_h :
        pyplot.plot([i,i])

pyplot.show()

print(c1*((z1+0.5*e1)**3 - (z1-0.5*e1)**3))
