print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP13 : COURBES ELLIPTIQUES                                                  #
# *************************************************************************** #
# *************************************************************************** #
""")

print("trabet clément")

# CONSIGNES
#
# Les seules lignes a modifier sont annoncee par "Code pour l'exercice"
# indique en commmentaire et son signalees
# Ne changez pas le nom des variables
#
# CONSEILS
#
# Ce modele vous sert a restituer votre travail. Il est deconseille d'ecrire
# une longue suite d'instruction et de debugger ensuite. Il vaut mieux tester
# le code que vous produisez ligne apres ligne, afficher les resultats et
# controler que les objets que vous definissez sont bien ceux que vous attendez.
#
# Vous devez verifier votre code en le testant, y compris par des exemples que
# vous aurez fabrique vous-meme.
#


reset()
print("""\
# ****************************************************************************
# ADDITION DANS UNE COURBE ELLIPTIQUE DONE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p = 61
Fp = FiniteField(p)
E = EllipticCurve(Fp,[1,0])
P = E.random_point()
Q = E.random_point()

# Code pour l'EXERCICE
def addition(P, Q):
    a = E.a4()
    xp, yp = Fp(P[0]), Fp(P[1])
    xq, yq = Fp(Q[0]), Fp(Q[1])
    lam=Fp(0)
    if P == O:
        return Q
    elif Q == 0:
        return P
    elif (xq == 0 and yq == 0):
        return P
    elif (xp == -xq and yp == yq):
        return O
    else:
        if (xq == xp and yq == yp):
            lam = (3 * xp**2+a) / (2 * yp)
        else:
            lam = (yq - yp) / (xq - xp)
        
        Rx = lam**2 - xp - xq
        Ry = lam * (xp - Rx) - yp
        R = E(Rx, Ry)

        return R


# # Affichage des resultats

addition(P,Q)


reset()
print("""\
# ****************************************************************************
# COURBE DE L'ANSSI DONE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice


# Code pour l'EXERCICE

ANSSI = "Agence nationale de la sécurité des systèmes d'information"
p = 109454571331697278617670725030735128146004546811402412653072203207726079563233
a = 109454571331697278617670725030735128145969349647868738157201323556196022393856
b = 107744541122042688792155207242782455150382764043089114141096634497567301547839
E = EllipticCurve(FiniteField(p),[a,b])

# # Affichage des resultats

print("ANSSI signifie :",ANSSI)
print("La courbe recommandée est")
print(E) 
print(p.is_prime())

print("En calculant 2**128 on remarque que la valeur est proche de racine de p")



reset()
print("""\
# ****************************************************************************
# COMPTAGE DE POINTS DONE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p = 2003
Fp = FiniteField(p)
a = 1929
b = 1178

while true:
    d=Fp.random_element()
    if not d.is_square():
        break


# Code pour l'EXERCICE

def comptage(a,b):
    Pol.<x> = PolynomialRing(Fp)
    points = set()
    for y in Fp :
        f = x^3 + a*x + b -y^2
        roots = f.roots()
        for xx in roots :
            points.add((xx[0],y))
    return len(points)+1 # il faut rajouter le point à l'infini

p = 11
Fp = FiniteField(p)    
    
frequence = [0 for _ in range(p+1+2*ceil(sqrt(p)))]
    
for aa in Fp :
    for bb in Fp :
        if -(4*aa^3 + 27*bb^2) != Fp(0) :
            n = EllipticCurve(Fp,[aa,bb]).cardinality()
            frequence[n]+=1

p = 2003
Fp = FiniteField(p)       
            
# # Affichage des resultats

print("Algorithme naif :", comptage(a,b))
print("Sage :", EllipticCurve(Fp,[a,b]).cardinality())
chart = bar_chart(frequence)
show(chart)


reset()
print("""\
# ****************************************************************************
# FACTORISATION ECM  DONE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

class FoundFactor(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

n = 2020

# Code pour l'EXERCICE

def division(x,y):
    try :
        quotient = x/y
        return quotient
    except ZeroDivisionError :
        raise FoundFactor(y)

def addition(P,Q):
    E = P.curve()
    Fp = E.base_ring()
    a = E.a4()
    b = E.a6()
    if P == 0 :
        return Q
    elif Q == 0 :
        return P
    try :
        xp, yp = P.xy()
    except ZeroDivisionError : # la troisième coordonnée projective du point n'est pas inversible
        factor = P[2]
        raise FoundFactor(factor)
    try :
        xq, yq = Q.xy()
    except ZeroDivisionError : # idem
        factor = Q[3]
        raise FoundFactor(factor)  
    if P != Q and P!=-Q :
        lamb = division((yq-yp),(xq-xp))
    elif P==Q :
        lamb = division((3*xp^2+a),2*yp)
    else : # on est dans le cas P = -Q
        return E(0)
    x = lamb^2 - xp - xq
    y = -yp + lamb*(xp-x)
    return E(x,y)

       
    
def multiplication(lamb,P): 
    E = P.curve()
    if lamb == 0 :
        return E(0)
    elif lamb%2 == 0:
        return multiplication(addition(P,P),lamb//2)
    else :
        return addition(multiplication(addition(P,P),lamb//2),P)

    

def ECM(n,B):
    Zn = Zmod(n)
    a, x0, y0 = Zn.random_element(), Zn.random_element(), Zn.random_element()
    b = y0^2 - x0^3 - a*x0
    g = gcd(4*a^3+27*b^2,n)
    if g > 1 and g < n :
        return g
    elif g == n :
        return None
    E = EllipticCurve(Zn,[a,b])
    A = E(x0,y0)
    for p in Primes(B) :
        e = 1
        while p^e <= B :
            e+=1
        try :
            A = multiplication(p^e,A)
        except FoundFactor as ex :
            factor = ex.value
            return gcd(factor,n)   
    return None


# # Affichage des resultats

print(ECM(n,15))

reset()
print("""\
# ****************************************************************************
# EXPONENTIATION RAPIDE ET ATTAQUE PAR CANAUX CACHES
# ****************************************************************************
""")

# NE PAS TRAITER



reset()
print("""\
# ****************************************************************************
# COURBE D'EDWARDS DONE
# ****************************************************************************
""")



# Code pour l'EXERCICE

reponse = "Il y'a différents avantages cryptographique a ce type de courbe, tout d'abord, la loi d'addition est la même pour tout les points de l'espace rendant l'algorithme résistant aux attaques par canaux auxiliaires, l'implémentation est donc moins complexe et les calculs plus rapides"

# # Affichage des resultats

print(reponse)

