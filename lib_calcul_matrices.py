
# coding: utf-8

#import scipy as sc
from sympy import *
def round_expr(expr, num_digits):
    return expr.xreplace({n : round(n, num_digits) for n in expr.atoms(Number)})
from pprint import pprint
from matplotlib import pyplot




def cbarre(E1,E2,v12,G12):
    '''calcul la matrice de rigidite a partir des grandeurs generalisees'''
    ###Determination de la matrice de rigidite a partir des grandeurs generalisees
    Q11 = E1/(1-E2/E1*v12**2)
    Q22 = E2/E1*Q11
    Q12 = v12*Q22
    Q66 = G12
    ###Matrice de rigidite dans la base d'orthropie
    cbarre = Matrix([Q11,Q12,0,
            Q12,Q22,0,
            0,0,Q66]).reshape(3,3)
    sbarre = Matrix([1/E1,-v12/E1,0,
            -v12/E1,1/E2,0,
            0,0,1/G12]).reshape(3,3)
    print("pour les valeurs de E1, E2, v12, G12 suivantes :   ")
    pprint([E1,E2,v12,G12])
    print()
    print("la matrice de rigidite en Gpa (10^9 Pa)dans la base d'orthropie = \n")
    pprint(round_expr(cbarre,2))
    print()
    print("la matrice de souplesse en Gpa (10^9 Pa)dans la base d'orthropie = \n")
    pprint(round_expr(sbarre,6))
    
    return cbarre





def ccouche(couche,alpha,cbarre):
    '''effectue un changement de base pour passer de la matrice de rigidite de la base d orthotropie vers les couches'''
    ###matrice de rigidite dans une nouvelle base
    c = cos(alpha)
    s = sin(alpha)
    cdbrigid = Matrix([c**4,s**4,2*c**2*s**2,4*c**2*s**2,
                    s**4,c**4,2*c**2*s**2,4*c**2*s**2,
                    c**2*s**2,c**2*s**2,c**4+s**4,-4*c**2*s**2,
                    c**2*s**2,c**2*s**2,-2*c**2*s**2,(c**2-s**2)**2,
                    c**3*s,-c*s**3,-c*s*(c**2-s**2),-2*c*s*(c**2-s**2),
                    c*s**3,-c**3*s,c*s*(c**2-s**2),2*c*s*(c**2-s**2)]).reshape(6,4)
    Q11 = cbarre[0]
    Q22 = cbarre[4]
    Q12 = cbarre[1]
    Q66 = cbarre[8]
    Q = Matrix([Q11,Q22,Q12,Q66]).reshape(4,1)
    C = cdbrigid*Q
    C1 = Matrix([C[0],C[2],C[4],C[2],C[1],C[5],C[4],C[5],C[3]]).reshape(3,3)
    print()
    print("Matrice de changement de base pour alpha = %s " % (couche))
    #pprint(round_expr(simplify(cdbrigid),2))
    pprint(simplify(cdbrigid))
    print()
    print("Matrice de rigidité en Gpa (10^9 Pa) pour alpha = %s " % (couche))
    pprint(round_expr(C1,2))
    return C1



def cdb_sigma(mat,alph) :
    '''effectue un changement de base pour passer des contraintes dans la couche vers les contraintes dans la base d orthotropie'''
    c = cos(alph)
    s = sin(alph)
    R = Matrix([
        [c**2, s**2 ,  2*c*s],
        [s**2 , c**2 , -2*c*s],
        [ -c*s, c*s, c**2-s**2]])
    print("Matrice de changement de base pour alpha = %s \n" % (alph))
    pprint(simplify(R))
    return R.evalf()*mat



def cont_max(mat_sig,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r ):
    '''calcul le coeficient de contrainte maximale a partir de la matrice de rigidite ramenee dans la base d orthotropie
    attention matrice de rigidite en Gpa et les contraintes a rupture en Mpa'''
    if mat_sig[0] > 0 : 
        F1_sig = SIG1_rt/mat_sig[0]
        print("Critere de contrainte maximale en traction : FSIG1 = %s "% (round(F1_sig,4)))
    else :
        F1_sig = SIG1_rc/mat_sig[0]
        print("Critere de contrainte maximale en compression : FSIG1 = %s "% (round(F1_sig,4)))

    if mat_sig[1] > 0 : 
        F2_sig = SIG2_rt/mat_sig[1]
        print("Critere de contrainte maximale en traction : FSIG2 = %s "% (round(F2_sig,4)))
    else :
        F2_sig = SIG2_rc/mat_sig[1]
        print("Critere de contrainte maximale en compression : FSIG2 = %s "% (round(F2_sig,4)))

    F6_sig = SIG6_r/mat_sig[2]
    print("Critere de contrainte maximale FSIG6  = %s"% (round(F6_sig,4)))
    return([F1_sig,F2_sig,F6_sig])


def tsai_hill(mat_sig,SIG1_rt,SIG1_rc,SIG2_rt,SIG2_rc,SIG6_r ):
    '''calcul le coeficient de contrainte tsai hill a partir de la matrice de rigidite ramenee dans la base d orthotropie
    attention matrice de rigidite en Mpa et les contraintes a rupture en Mpa'''
    SIG1 = mat_sig[0]
    SIG2 = mat_sig[1]
    SIG6 = mat_sig[2]

    if SIG1 > 0 :
        F11 = SIG1**2/SIG1_rt**2

    else :
        F11 = SIG1**2/SIG1_rc**2

    print("F11 = %s"% (round(F11,8)))

    if SIG2 > 0 :
        F22 = SIG2**2/SIG2_rt**2

    else :
        F22 = SIG2**2/SIG2_rc**2

    print("F22 = %s"% (round(F22,8)))
    
    if SIG1 > 0 :
        F12 = SIG1*SIG2/SIG1_rt**2

    else :
        F12 = SIG1*SIG2/SIG1_rc**2
    print("F12 = %s"% (round(F12,8)))
   
    F66 = SIG6**2/SIG6_r**2

    print("F66 = %s"% (round(F66,8)))

    F_th = 1/sqrt(F11+F22+F12+F66)
    pprint("Critere de contrainte Tsai hill %s"% (round(F_th,4)))
   
    
