import numpy as np

LondonEyeSpokeRadius = 0.055
LondonEyeSpokeArea = np.pi*LondonEyeSpokeRadius**2
EnterpriseSpokeRadius = 0.055 * 0.135
EnterpriseSpokeArea = np.pi*EnterpriseSpokeRadius**2
YieldStrength = 1180*10**6



def tensions_function_mass(m, omega, r, g, T_0 ):
    masses = [m, m*2, m*4, m*6, m*8]
    Tensions = []
    Stresses = []
    FOSs = []
    for mass in masses:
        T = T_0 + mass*omega**2*r - mass*g*np.sin(-np.pi/2)
        Tensions.append(T)
        Stresses.append(T/LondonEyeSpokeArea)
        FOSs.append(YieldStrength/(T/LondonEyeSpokeArea))
    print(f"Values of mass: {masses}")
    print(f"London Eye Tensions, MASS [m, m*2, m*4, m*6, m*8]: {Tensions}")
    print(f"London Eye Stresses, MASS [m, m*2, m*4, m*6, m*8]: {Stresses}")
    print(f"London Eye FOSs, MASS [m, m*2, m*4, m*6, m*8]: {FOSs}")
    return Tensions
   


def tensions_function_radius(m, omega, r, g, T_0 ):
    radiuses = [r, r*2, r*10, r*100, r*1000]
    Tensions = []
    Stresses = []
    FOSs = []
    for radius in radiuses:
        T = T_0 + m*omega**2*radius - m*g*np.sin(-np.pi/2)
        Tensions.append(T)
        Stresses.append(T/LondonEyeSpokeArea)
        FOSs.append(YieldStrength/(T/LondonEyeSpokeArea))
    print(f"Values of radius: {radiuses}")
    print(f"London Eye Tensions, RADIUS [r, r*2, r*10, r*100, r*1000]: {Tensions}")
    print(f"London Eye Stresses, RADIUS [r, r*2, r*10, r*100, r*1000]: {Stresses}")
    print(f"London Eye FOSs, RADIUS [r, r*2, r*10, r*100, r*1000]: {FOSs}")
    return Tensions
   
                            


def tensions_function_omega(m, omega, r, g, T_0 ):
    omegas = [omega, omega*5, omega*10, omega*50, omega*100]
    Tensions = []
    Stresses = []
    FOSs = []
    for angle in omegas:
        T = T_0 + m*angle**2*r - m*g*np.sin(-np.pi/2)
        Tensions.append(T)
        Stresses.append(T/LondonEyeSpokeArea)
        FOSs.append(YieldStrength/(T/LondonEyeSpokeArea))
    print(f"Values of Omega: {omegas}")
    print(f"London Eye Tensions [omega, omega*5, omega*10, omega*50, omega*100]: {Tensions}")
    print(f"London Eye Stresses [omega, omega*5, omega*10, omega*50, omega*100]: {Stresses}")
    print(f"London Eye FOSs [omega, omega*5, omega*10, omega*50, omega*100]: {FOSs}")
    return Tensions





def Etensions_function_mass(m, omega, r, g, T_0 ):
    masses = [m, m*1.25, m*1.5, m*1.75, m*2]
    Tensions = []
    Stresses = []
    FOSs = []
    for mass in masses:
        T = T_0 + mass*omega**2*r - mass*g*np.sin(-np.pi/2)
        Tensions.append(T)
        Stresses.append(T/EnterpriseSpokeArea)
        FOSs.append(YieldStrength/(T/EnterpriseSpokeArea))
    print(f"Values of mass: {masses}")
    print(f"Enterprise Tensions, MASS [m, m*1.75, m*2.5, m*3.25, m*4]: {Tensions}")
    print(f"Enterprise Stresses, MASS [m, m*1.75, m*2.5, m*3.25, m*4]: {Stresses}")
    print(f"Enterprisee FOSs, MASS [m, m*1.75, m*2.5, m*3.25, m*4]: {FOSs}")
    return Tensions
   


def Etensions_function_radius(m, omega, r, g, T_0 ):
    radiuses = [r, r*1.5, r*2, r*2.5, r*3]
    Tensions = []
    Stresses = []
    FOSs = []
    for radius in radiuses:
        T = T_0 + m*omega**2*radius - m*g*np.sin(-np.pi/2)
        Tensions.append(T)
        Stresses.append(T/EnterpriseSpokeArea)
        FOSs.append(YieldStrength/(T/EnterpriseSpokeArea))
    print(f"Values of radius: {radiuses}")
    print(f"Enterprise Tensions, RADIUS [r, r*2, r*3, r*4, r*5]: {Tensions}")
    print(f"Enterprise Stresses, RADIUS [r, r*2, r*3, r*4, r*5]: {Stresses}")
    print(f"Enterprise FOSs, RADIUS [r, r*2, r*3, r*4, r*5]: {FOSs}")
    return Tensions
   
                            


def Etensions_function_omega(m, omega, r, g, T_0 ):
    omegas = [omega, omega*1.25, omega*1.5, omega*1.75, omega*2]
    Tensions = []
    Stresses = []
    FOSs = []
    for angle in omegas:
        T = T_0 + m*angle**2*r - m*g*np.sin(-np.pi/2)
        Tensions.append(T)
        Stresses.append(T/EnterpriseSpokeArea)
        FOSs.append(YieldStrength/(T/EnterpriseSpokeArea))
    print(f"Values of Omega: {omegas}")
    print(f"Enterprise Tensions [omega, omega*1.25, omega*1.5, omega*1.75, omega*2]: {Tensions}")
    print(f"Enterprise Stresses [omega, omega*1.25, omega*1.5, omega*1.75, omega*2]: {Stresses}")
    print(f"Enterprise FOSs [omega, omega*1.25, omega*1.5, omega*1.75, omega*2]: {FOSs}")
    return Tensions


LEPre_Tension = 26525*(9.81 - 0.00349**2*60)
EnterprisePre_Tension = 1000

# Print Tension, Stress & FOS for parameter values
tensions_function_mass(26525, 0.00349, 60, 9.81, LEPre_Tension)
print()
tensions_function_radius(26525, 0.00349, 60, 9.81, LEPre_Tension)
print()
tensions_function_omega(26525, 0.00349, 60, 9.81, LEPre_Tension)
print()
print()
Etensions_function_mass(2000, 1.414, 8.1, 9.81, EnterprisePre_Tension)
print()
Etensions_function_radius(2000, 1.414, 8.1, 9.81, EnterprisePre_Tension)
print()
Etensions_function_omega(2000, 1.414, 8.1, 9.81, EnterprisePre_Tension)