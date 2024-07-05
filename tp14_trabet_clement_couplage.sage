print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP14 : LOG DISCRET ET COUPLAGES                                             #
# *************************************************************************** #
# *************************************************************************** #
""")

print("TRABET Clément")
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
# PAS DE BEBE, PAS DE GEANT DONE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p1 = 1823
Fp1 = FiniteField(p1)
b1 = Fp1(3)
x1 = Fp1(693)

p2 = 239
Fp2 = FiniteField(p2)
b2 = Fp2(2)
x2 = Fp2(15)


# Code pour l'EXERCICE

def Shanks(x,b):
    Fp = x.parent()
    p = Fp.cardinality()
    s = ceil(sqrt(p-1))
    T = dict()
    for j in range(s) :
        beta = x*b^(-j)
        T[beta] = j
    i = 0
    gamma = Fp(1)
    bprime = b^s
    while gamma not in T.keys() :
        i+=1
        gamma = gamma * bprime
    j = T[gamma]
    return i*s + j


# # Affichage des resultats

print("Question 2 :", Shanks(x1,b1))
print("Question 3 :", Shanks(x2,b2))




reset()
print("""\
# ****************************************************************************
# RHO DE POLLARD DONE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p= 281
Fp = FiniteField(p)
x1 = Fp(263) 
b1 = Fp(239)
x2 = Fp(165)
b2 = Fp(127)
x3 = Fp(210)
b3 = Fp(199)


# Code pour l'EXERCICE
def rho(g,b):
    Fp = g.parent()
    p = Fp.cardinality()
    Zn = Zmod(p-1)
    partition = lambda x : hash(x)%3
    def phi(w,alpha,beta) :
        if partition(w) == 0 :
            return (g*w, alpha, beta+1)
        elif partition(w) == 1 :
            return (w^2, 2*alpha, 2*beta)
        elif partition(w) == 2 :
            return (b*w, alpha+1, beta)
    alpha, beta = randint(1,p), randint(1,p)
    x, ax, bx = phi(b^alpha*g^beta, alpha, beta)
    y, ay, by = phi(x,ax,bx)
    while x != y :
        x, ax, bx = phi(x ,ax ,bx)
        y, ay, by = phi(y, ay, by)
        y, ay, by = phi(y, ay, by)
    try :
        return Zn((ax-ay))/Zn((by-bx))
    except ZeroDivisionError :
        return rho(g,b)


# # Affichage des resultats

print("Le log de x=",x1,"en base",b1,"vaut",rho(x1,b1),".")
print("Le log de x=",x2,"en base",b2,"vaut",rho(x2,b2),".")
print("Le log de x=",x3,"en base",b3,"vaut",rho(x3,b3),".")





reset()
print("""\
# ****************************************************************************
# COUPLAGE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p = 61
Fp = FiniteField(p)
E = EllipticCurve(Fp,[11,0])
print("groupe de E=", E.abelian_group()) # pour verifier
S = E(24,34)
T = E(5,27)
r = 10
print("Verification de la r-torsion : r*S =", r*S, "et r*T =", r*T)


# Code pour l'EXERCICE


def myLine(P1,P2,S):
    E=P1.curve()
    x1=P1[0]; y1=P1[1]; z1=P1[2]
    x2=P2[0]; y2=P2[1]; z2=P2[2]
    xS=S[0]; yS=S[1]; zS=S[2]
    a = E.a4()
    b = E.a6()
    
    if P1 == P2 and y1 != 0:
        lamb = (3 * x1^2 + a) / (2 * y1)
        return yS - lamb*(xS - x1) - y1
    
    if P1 == P2 and y1 == 0:
        return z1*xS + x1*zS
    
    u = y1*z2 - y2*z1
    v = z1*x2 - z2*x1
    w = x1*y2 - x2*y1
    
    return u*xS + v*yS + w*zS
    
def myH(P1,P2,S):
    return myLine(P1,P2,S)/myLine(P1+P2,-P1-P2,S)

def myMiller(r,S,P):
    R = S
    f = 1
    b = r.bits()
    for r in b[-2::-1] :
        f = f^2 * myH(R,R,P)
        R = 2*R
        if r == 1 :
            f = f*myH(R,S,P)
            R = R + S
    return f

def myTatePairing(S,T,r):
    E = S.curve()
    while True :
        Q = E.random_element()
        if r*Q != E(0) :
            continue
        try :
            return myMiller(r,S,T+Q) / myMiller(r,S,Q)
        except :
            continue

# # Affichage des resultats

print("Calcul du couplage", myTatePairing(S,T,r))



reset()
print("""\
# ****************************************************************************
# ATTAQUE M.O.V.
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p = 2199023255579
Fp = FiniteField(p)
E = EllipticCurve(Fp,[1,0])
P = E(1435967701832 , 123951463462)
Q = E(1129476910351 , 1383670460733)

# Code pour l'EXERCICE

j = 0 # j-invariant a faire calculer par une fonction de SageMath
rep2 = "commentaire"
t = 2 # Ecrire le code pour calculer cette valeur
q = p^t
Fq.<alpha> = FiniteField(q)
EE = EllipticCurve(Fq,[1,0])
PP = EE(1435967701832 , 123951463462)
QQ = EE(1129476910351 , 1383670460733)
SS = EE(0) # point a calculer vous-meme
zeta1 = 1
zeta2 = 1
lambd = 1


# # Affichage des resultats

print("p premier ?",p.is_prime())
print("j-invariant de E :",j)
print("p mod 4 =", mod(p,4))
print(rep2)
print("Cardinal de E(Fp) :",E.cardinality(),"=",E.cardinality().factor())
print("Ordre de P :",P.order())
print("Cardinal de E(Fq) :",EE.cardinality(),"=",EE.cardinality().factor())
print("Point S :",SS)
print("On calcule zeta1 =",zeta1,", zeta2 =",zeta2,", lambda =",lambd,".")




reset()
print("""\
# ****************************************************************************
# CALCUL D'INDICE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

p = 439
Fp = FiniteField(p)
g = Fp(237)
b = Fp(136)
y = 11

def facto_friable(x,primes_y):
    fact = factor(x)
    expo = {}
    for pp,exp in fact :
        expo[pp] = exp
    for pp in primes_y :
        expo.setdefault(pp,0)
    vi = []
    for pp in primes_y :
        vi.append(expo[pp])
    return vi
    
def est_friable(x,primes_y):
    return all([pp in primes_y for (pp,_) in factor(x)])

def LogIndice(g,b,y):
    Fp = g.parent()
    p = Fp.cardinality()
    Zn = Zmod(p-1)
    primes = Primes()
    q = primes.first()
    primes_y = set()
    while q <= y :
        primes_y.add(q)
        q = primes.next(q)
    
    k = len(primes_y)
    M = (ZZ^k)/((p-1)*ZZ^k)
    alphas = []
    v = []
    i = 1
    while i <= 4*k :
        alpha = randint(0,p-2)
        gamma = int(b^alpha)
        if est_friable(gamma,primes_y) :
            vi = facto_friable(gamma,primes_y)
            v.append(vi)
            alphas.append(alpha)
            i+=1
    while M.submodule(v) != M :
        alphas = []
        v = []
        i = 1
        while i <= 4*k :
            alpha = randint(0,p-2)
            gamma = int(b^alpha)
            if est_friable(gamma,primes_y) :
                vi = facto_friable(gamma,primes_y)
                v.append(vi)
                alphas.append(alpha)
                i+=1
        
    M = matrix(Zn,v)
    alphas = vector(Zn,alphas)
    logs = M.solve_right(alphas)
    beta = randint(0,p-2)
    while not est_friable(int(b^beta*g),primes_y) :
        beta = randint(0,p-2)
    facto = facto_friable(int(b^beta*g),primes_y)
    return -beta + sum([logs[i]*facto[i] for i in range(k)])

# # Affichage des resultats

print("Le log de g=",g,"en base",b,"vaut",LogIndice(g,b,y),".")

