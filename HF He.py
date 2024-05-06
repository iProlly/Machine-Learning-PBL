import sympy as sym
import math as m
import re
import time



c1, c2, c3, c4, c5, c6, r, r1, r2 = sym.symbols("c1 c2 c3 c4 c5 c6 r r1 r2")
e = 2.71828

cut = [2,2,8,8,18,18,32]
s_con = []
spin=[]
orb_0=[0,0,0,0,0,0,0,0]
orb_1=[0,0,0,0,0,0,0,0]
orb_2=[0,0,0,0,0,0,0,0]
orb_3=[0,0,0,0,0,0,0,0]
proton = 2
f = proton

#CUTTING
for k in range(len(cut)):
    if proton - cut[k] <= 0:
        s_con.append(proton)
        break
    proton = proton - cut[k]
    s_con.append(cut[k])

#assign electron to orbital and assign spin
for j in range(len(s_con)):
    for p in range(j+2-m.ceil((j+1)/2), j+2):
        if j + 1 - p == 0:
            if f - 2 <= 0:
                globals()["orb_"+str(j+1-p)][p-1] += f
                for q in range(f):
                    if q+1 <= 1:
                        spin.append("+")
                    else:
                        spin.append("-")
                break
            globals()["orb_"+str(j+1-p)][p-1] += 2
            for q in range(2):
                    if q+1 <= 1:
                        spin.append("+")
                    else:
                        spin.append("-")
            f-=2
        if j + 1 - p == 1:
            if f - 6 <= 0:
                globals()["orb_"+str(j+1-p)][p-1] += f
                for q in range(f):
                    if q+1 <= 3:
                        spin.append("+")
                    else:
                        spin.append("-")
                break
            globals()["orb_"+str(j+1-p)][p-1] += 6
            for q in range(6):
                    if q+1 <= 3:
                        spin.append("+")
                    else:
                        spin.append("-")
            f-=6
        if j + 1 - p == 2:
            if f - 10 <= 0:
                globals()["orb_"+str(j+1-p)][p-1] += f
                for q in range(f):
                    if q+1 <= 5:
                        spin.append("+")
                    else:
                        spin.append("-")
                break
            globals()["orb_"+str(j+1-p)][p-1] += 10
            for q in range(10):
                    if q+1 <= 5:
                        spin.append("+")
                    else:
                        spin.append("-")
            f-=10
        if j + 1 - p == 3:
            if f - 14 <= 0:
                globals()["orb_"+str(j+1-p)][p-1] += f
                for q in range(f):
                    if q+1 <= 7:
                        spin.append("+")
                    else:
                        spin.append("-")
                break
            globals()["orb_"+str(j+1-p)][p-1] += 14
            for q in range(14):
                    if q+1 <= 7:
                        spin.append("+")
                    else:
                        spin.append("-")
            f-=14

#Alpha and Coefficient Database
alpha2 = [2.432879, 0.433051]
coeff2 = [0.430128, 0.678914]
alpha3 = [6.362421, 1.158923, 0.313650]
coeff3 = [0.154329, 0.535328, 0.444635]
alpha4 = [14.899830, 2.726485, 0.757447, 0.251390]
coeff4 = [0.0567524, 0.260141, 0.532846, 0.291625]
alpha5 = [32.290030, 5.917063, 1.652678, 0.564287, 0.212644]
coeff5 = [0.0221406, 0.113541, 0.331816, 0.482570, 0.193572]
alpha6 = [65.984568, 12.098198, 3.384640, 1.162715, 0.451516, 0.185959]
coeff6 = [0.00916360, 0.0493615, 0.168538, 0.370563, 0.416492, 0.130334]

#Construct wavefunction of each electron
psi = [0, 0]
n = input("put the number of primitive GTO: ")
start_time = time.time()
for i in range(2):
    for j in range(int(n)):
        psi[i] = psi[i] + globals()["coeff" + str(n)][j]*(2*globals()["alpha" + str(n)][j]/3.14159)**0.75*e**(-globals()["alpha" + str(n)][j]*globals()["r" + str(i+1)]**2)

#Evaluation of Gaussian integral with infinite upper limit
def gaussianintegral(z, dummy):
    #Split Algebraic Terms
    z = sym.expand(z)
    z_split = re.split('(\+ | -)', str(z))
    z_new = []
    z_new.append(sym.sympify(z_split[0]))
    for i in range(1, len(z_split), 2):
        z_new.append(sym.sympify(z_split[i]+z_split[i+1]))
    
    #Determine Parameters -> Use the formula
    sum = 0
    for j in range(len(z_new)):
        if (str(z_new[j])[str(z_new[j]).rfind('/')-2]=='r'):
            power = 1
        elif (str(z_new[j])[str(z_new[j]).rfind('/')-2]!='*'):
            power = 0
        else:
            power = sym.sympify(str(z_new[j])[str(z_new[j]).rfind('/')-1])
        coeff = (z_new[j]/sym.sympify(dummy)**power).subs(dummy, 0)
        expo = m.log(((z_new[j]/sym.sympify(dummy)**power).subs(dummy, 0))/((z_new[j]/sym.sympify(dummy)**power).subs(dummy, 1)))
        
        sum += coeff*m.gamma((power+1)/2)/(2*expo**((power+1)/2))
    return sum

#Evaluation of Integral contain erf
def sperf2integral(z, dummy):
    #Split Algebraic Terms
    z = sym.expand(z)
    z_split = re.split('(\+ | -)', str(z))
    z_new = []
    z_new.append(sym.sympify(z_split[0]))
    for i in range(1, len(z_split), 2):
        z_new.append(sym.sympify(z_split[i]+z_split[i+1]))
    
    sum = 0
    for j in range(len(z_new)):
        #Check if there is erf or not
        if "erf" in str(z_new[j]):
            #Determine Parameters -> Use the formula
            z_new[j] = sym.expand(z_new[j]/dummy)
            z_new_split = list(filter(None, re.split("[*()/]", str(z_new[j]))))
            sum += float(z_new_split[0])*float(z_new_split[2])/(2*float(z_new_split[5])*(float(z_new_split[2])**2+float(z_new_split[5]))**0.5)
        else:
            #Determine Parameters -> Use the formula
            if (str(z_new[j])[str(z_new[j]).rfind('/')-2]!='*'):
                power = 0
            else:
                power = sym.sympify(str(z_new[j])[str(z_new[j]).rfind('/')-1])
            coeff = (z_new[j]/sym.sympify(dummy)**power).subs(dummy, 0)
            expo = m.log(((z_new[j]/sym.sympify(dummy)**power).subs(dummy, 0))/((z_new[j]/sym.sympify(dummy)**power).subs(dummy, 1)))   
            if (int(power) == 2):
                sum += coeff*(m.pi)**0.5/(4*expo**1.5)
    return sum

#Evaluation of Gaussian integral with finite upper limit
def Gaussianintegralfin(z, dummy, bound):
    #Split Algebraic Terms
    z = sym.expand(z)
    z_split = re.split('(\+ | -)', str(z))
    z_new = []
    z_new.append(sym.sympify(z_split[0]))
    for i in range(1, len(z_split), 2):
        z_new.append(sym.sympify(z_split[i]+z_split[i+1]))
    
    #Determine Parameters -> Use the formula
    sum = 0
    for j in range(len(z_new)): 
        if (str(z_new[j])[str(z_new[j]).rfind('/')-2]!='*'):
            power = 0
        else:
            power = sym.sympify(str(z_new[j])[str(z_new[j]).rfind('/')-1])
        coeff = (z_new[j]/sym.sympify(dummy)**power).subs(dummy, 0)
        expo = m.log(((z_new[j]/sym.sympify(dummy)**power).subs(dummy, 0))/((z_new[j]/sym.sympify(dummy)**power).subs(dummy, 1)))   
        if (int(power) == 2):
            sum += -coeff/(2*expo)*(bound*e**(-expo*bound**2)-m.pi**0.5*sym.erf(expo**0.5*bound)/(2*expo**0.5)) 
    return sum

#Evaluation of Kinetic energy
def t(psi, i):
    return (gaussianintegral(globals()["r" + str(i)]**2*psi[i-1]*(-1/(2*globals()["r" + str(i)])*sym.diff(globals()["r" + str(i)]*psi[i-1], globals()["r" + str(i)], 2)), globals()["r" + str(i)]))/(gaussianintegral(globals()["r" + str(i)]**2*psi[i-1]**2, globals()["r" + str(i)]))

#Evaluation of electron-nucleus interaction
def v(psi, i):
    return (gaussianintegral(globals()["r" + str(i)]**2*psi[i-1]*(-2*psi[i-1]/globals()["r" + str(i)]), globals()["r" + str(i)]))/(gaussianintegral(globals()["r" + str(i)]**2*psi[i-1]**2, globals()["r" + str(i)]))

#Evaluation of electron-electron interaction
def j(psi, i):
    integral = 0
    for k in range(1, 3):
        integral1 = 0
        integral2 = 0
        if (k != i):
            integral1 += 16*3.14159**2*(sperf2integral(psi[i-1]**2*globals()["r" + str(i)]*Gaussianintegralfin(psi[k-1]**2*globals()["r" + str(k)]**2, globals()["r" + str(k)], globals()["r" + str(i)]), globals()["r" + str(i)]))
            integral2 += 16*3.14159**2*(sperf2integral(psi[k-1]**2*globals()["r" + str(k)]*Gaussianintegralfin(psi[i-1]**2*globals()["r" + str(i)]**2, globals()["r" + str(i)], globals()["r" + str(k)]), globals()["r" + str(k)]))
        integral += 0.5*(integral1+integral2)
    return integral

#Evaluation of exchange energy
def k(psi, i):
    integral = 0
    for k in range(1, 3):
        integral1 = 0
        integral2 = 0
        if (k != i):
            if spin[k-1] == spin[i-1]:
                integral1 += 16*3.14159**2*(sperf2integral(psi[i-1]*psi[k-1].subs(globals()["r" + str(k)], globals()["r" + str(i)])*globals()["r" + str(i)]*Gaussianintegralfin(psi[k-1]*psi[i-1].subs(globals()["r" + str(i)], globals()["r" + str(k)])*globals()["r" + str(k)]**2, globals()["r" + str(k)], globals()["r" + str(i)]), globals()["r" + str(i)]))
                integral2 += 16*3.14159**2*(sperf2integral(psi[k-1]*psi[i-1].subs(globals()["r" + str(i)], globals()["r" + str(k)])*globals()["r" + str(k)]*Gaussianintegralfin(psi[i-1]*psi[k-1].subs(globals()["r" + str(k)], globals()["r" + str(i)])*globals()["r" + str(i)]**2, globals()["r" + str(i)], globals()["r" + str(k)]), globals()["r" + str(k)]))
            else:
                integral1 += 0
                integral2 += 0
        integral += 0.5*(integral1+integral2)
    return integral

#Evaluation of Overall energy
def E(psi, i):
    return t(psi, i) + v(psi, i) + j(psi, i) - k(psi, i)

print("Helium atom Energy from STO-"+str(n)+"G using Hartree-Fock Method= ", E(psi,1)+E(psi,2))
print("Time used: %s seconds" % (time.time() - start_time))
