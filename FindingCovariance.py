'''
 This routine is part of the publication 'The Impact of
 Magnetospheric Currents on Jovian Magnetic Field Reconstruction' 
 (Loncar & Jackson, 2026) and developed to find the covariance
 between magentic field measurements as a result of stochastic 
 magnetospheric currents. 
 
 Copyright (C) 2026 Marco Loncar

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.

########################################################################

 The structure of this code mirrors that of the mathematical 
 procedure used in the associated paper (with relevant equations 
 referenced in parentheses) and so is not necessarily designed 
 for efficiency but rather for ease of understanding.

 Covariances are found as follows:
     - Input measurment positions
     - Check for r < r' as required for (14) to be valid
     - Find angular separation and (numerically) calculate Legendre
        polynomials and their derivatives
     - Using the equations of (B5), calculate the covariance of 
        desired field components
     - Derivatives of Upsilon are explicitly found through use of 
        the chain rule and the truncated derivatives of (14)
'''

''' =======================================
                  LIBRARIES
    ======================================= '''

import numpy as np
import scipy as sp

''' =======================================
                  PARAMETERS
    ======================================= '''

# Constants
#   Planetary radius (currently taken to be Jupiter)
R_P = 71492000
#   Free space permeability
mu0 = sp.constants.mu_0

# Current correlation params (as appearing in the associated paper)
#   Radial correlation width, alpha (units of R_P)
alpha = 0.01
#   Angular correlation control, gamma
gamma = 0.9
#   Inner and outer magnetosphere radii (units of R_P)
RI = 1; RO = 6

# Sum truncation degree
Ntrunc = 200

# Mu cut-off near +-1
#   [at small angular separations, this ensures computational
#    error does lead to incorrect evaluation of mu]
CutOff = 1e-5

''' =======================================
                  FUNCTIONS
    ======================================= '''

# Calculate P_n(x), P'_n(x), P''_n(x) up to Nmax
def LegP_finder(x):
    # Initialise arrays for Legendre polynomials and their derivatives
    LegP = np.zeros(Ntrunc + 2)
    dLegP = np.zeros(Ntrunc)
    ddLegP = np.zeros(Ntrunc)
    
    # Calcuate the first terms of P_n and P'_n directly
    LegP[0] = 1; LegP[1] = x; dLegP[1] = 1
    
    # Expressions for special cases x = +-1 (analytically found)
    if (x == 1):
        for n in range(Ntrunc):
            LegP = np.ones(Ntrunc + 2)
            dLegP[n] = 0.5 * n * (n + 1)
            ddLegP[n] = 0.125 * (n - 1) * n * (n + 1) * (n + 2)
    elif (x == -1):
        for n in range(Ntrunc):
            LegP[n] = (-1)**n
            dLegP[n] = 0.5 * n * (n + 1) * (-1)**(n + 1)
            ddLegP[n] = 0.125 * (n - 1) * n * (n + 1) * (n + 2) * (-1)**n
    # Recursive forms of LegP elements for x != +-1
    else:
        # Fill the rest of the P_n(x) array
        for n in range(Ntrunc):
            # Use the recurrence relation for P_n+2(x)
            LegP[n+2] = ((2 * n + 3) * x * LegP[n+1] - (n + 1) * LegP[n]) / (n + 2)
        # Use this array to find the derivatives of P_n(x)
        for n in range(2, Ntrunc):
            # Using the recurrence relations for derivatives
            dLegP[n] = (n + 1) * (LegP[n+1] - x * LegP[n]) / (x**2 - 1)
            ddLegP[n] = (n + 1) * ((n + 2) * LegP[n+2] - (2 * n + 5) * x * LegP[n+1] + ((n + 2) * (x**2) + 1) * LegP[n]) / ((x**2 - 1)**2)
    
    # Output results
    return LegP, dLegP, ddLegP

#############################################################################

''' 
    To reduce (some of) the coming clutter, mu and its 
    derivatives are explicitly defined here 
        - mu is a measure of separation between 
            two measurements' angular positions:
                > (theta, phi) = (t, p)
                > (theta', phi') = (T, P)
        - mu = cos(angular separation)
        - -1 <= mu <= 1
'''

def dmu_dt(t, T, p, P):
    result = np.cos(t) * np.sin(T) * np.cos(p - P) - np.sin(t) * np.cos(T)
    return result
def dmu_dT(t, T, p, P):
    result = np.sin(t) * np.cos(T) * np.cos(p - P) - np.cos(t) * np.sin(T)
    return result
def dmu_dp(t, T, p, P):
    result = - np.sin(t) * np.sin(T) * np.sin(p - P)
    return result
def dmu_dP(t, T, p, P):
    result = np.sin(t) * np.sin(T) * np.sin(p - P)
    return result
def d2mu_dtdT(t, T, p, P):
    result = np.cos(t) * np.cos(T) * np.cos(p - P) + np.sin(t) * np.sin(T)
    return result
def d2mu_dtdP(t, T, p, P):
    result = np.cos(t) * np.sin(T) * np.sin(p - P)
    return result
def d2mu_dpdT(t, T, p, P):
    result = - np.sin(t) * np.cos(T) * np.sin(p - P)
    return result
def d2mu_dpdP(t, T, p, P):
    result = np.sin(t) * np.sin(T) * np.cos(p - P)
    return result

#############################################################################

'''
    Upsilon derivatives
        - Derivatives of the sum given in (14)
        - Separated into 'Inner', 'Outer' and 'Mid' contributions
            [corresponding to the three sums of (14)]
        - Summed up to truncation degree Ntrunc
        - The three contributions are appropriately scaled and 
            combined to give the full result
'''

# Calculate double radial (r,R) derivative
def d2Y_drdR(r, R, Leg):
    # Intialising contributions from boundaries and centre
    InnerCont = 0; OuterCont = 0; MidCont = 0
    for n in range(Ntrunc):
        InnerCont += (n + 1)**2 * (gamma * (RI**2) / (r * R))**(n + 1) * Leg[n] / ((2 * n + 1) * (2 * n + 5) * r * R)
        OuterCont += n**2 * ((gamma * r * R) / (RO**2))**n * Leg[n] / ((2 * n + 1) * (2 * n - 3) * r * R)
        MidCont += (gamma * r / R)**n * Leg[n] * ((n * (n - 3) / (2 * n - 3)) - ((n + 1) * (n + 4) / (2 * n + 5)) * ((r / R)**4))
    # Final scaling
    InnerCont *= -48 * np.pi * (RI**3) / (gamma * (alpha**2))
    OuterCont *= -48 * np.pi * (RO**3) / (alpha**2)
    MidCont *= -12 * np.pi * (R**2) / (r * alpha**2)
    
    # Add results together
    result = InnerCont + OuterCont + MidCont
    return result

# Calculate mixed radial-angular (r,mu) derivative
def d2Y_drdmu(r, R, dLeg):
    # Intialising contributions from boundaries and centre
    InnerCont = 0; OuterCont = 0; MidCont = 0
    for n in range(Ntrunc):
        InnerCont += - (n + 1) * (gamma * (RI**2) / (r * R))**(n + 1) * dLeg[n] / ((2 * n + 1) * (2 * n + 5) * r)
        OuterCont += n * ((gamma * r * R) / (RO**2))**n * dLeg[n] / ((2 * n + 1) * (2 * n - 3) * r)
        MidCont += (gamma * r / R)**n * dLeg[n] * ((n / (2 * n - 3)) - ((n + 4) / (2 * n + 5)) * ((r / R)**4))
    # Final scaling
    InnerCont *= -48 * np.pi * (RI**3) / (gamma * (alpha**2))
    OuterCont *= -48 * np.pi * (RO**3) / (alpha**2)
    MidCont *= 12 * np.pi * (R**3) / (r * alpha**2)
    
    # Add results together
    result = InnerCont + OuterCont + MidCont
    return result

# Calculate mixed radial-angular (r',mu) derivative
def d2Y_dRdmu(r, R, dLeg):
    # Intialising contributions from boundaries and centre
    InnerCont = 0; OuterCont = 0; MidCont = 0
    for n in range(Ntrunc):
        InnerCont += - (n + 1) * (gamma * (RI**2) / (r * R))**(n + 1) * dLeg[n] / ((2 * n + 1) * (2 * n + 5) * R)
        OuterCont += n * ((gamma * r * R) / (RO**2))**n * dLeg[n] / ((2 * n + 1) * (2 * n - 3) * R)
        MidCont += (gamma * r / R)**n * dLeg[n] * (((n - 3) / (2 * n - 3)) - ((n + 1) / (2 * n + 5)) * ((r / R)**4))
    # Final scaling
    InnerCont *= -48 * np.pi * (RI**3) / (gamma * (alpha**2))
    OuterCont *= -48 * np.pi * (RO**3) / (alpha**2)
    MidCont *= -12 * np.pi * (R**2) / (alpha**2)
    
    # Add results together
    result = InnerCont + OuterCont + MidCont
    return result

# Calculate single angular (mu) derivative
def dY_dmu(r, R, dLeg):
    # Intialising contributions from boundaries and centre
    InnerCont = 0; OuterCont = 0; MidCont = 0
    for n in range(Ntrunc):
        InnerCont += (gamma * (RI**2) / (r * R))**(n + 1) * dLeg[n] / ((2 * n + 1) * (2 * n + 5))
        OuterCont += ((gamma * r * R) / (RO**2))**n * dLeg[n] / ((2 * n + 1) * (2 * n - 3))
        MidCont += (gamma * r / R)**n * dLeg[n] * ((1 / (2 * n - 3)) - (1 / (2 * n + 5)) * ((r / R)**4))
    # Final scaling
    InnerCont *= -48 * np.pi * (RI**3) / (gamma * (alpha**2))
    OuterCont *= -48 * np.pi * (RO**3) / (alpha**2)
    MidCont *= 12 * np.pi * (R**3) / (alpha**2)
    
    # Add results together
    result = InnerCont + OuterCont + MidCont
    return result
    
# Calculate second angular (mu,mu) derivative
def d2Y_ddmu(r, R, ddLeg):
    # Intialising contributions from boundaries and centre
    InnerCont = 0; OuterCont = 0; MidCont = 0
    for n in range(Ntrunc):
        InnerCont += (gamma * (RI**2) / (r * R))**(n + 1) * ddLeg[n] / ((2 * n + 1) * (2 * n + 5))
        OuterCont += ((gamma * r * R) / (RO**2))**n * ddLeg[n] / ((2 * n + 1) * (2 * n - 3))
        MidCont += (gamma * r / R)**n * ddLeg[n] * ((1 / (2 * n - 3)) - (1 / (2 * n + 5)) * ((r / R)**4))
    # Final scaling
    InnerCont *= -48 * np.pi * (RI**3) / (gamma * (alpha**2))
    OuterCont *= -48 * np.pi * (RO**3) / (alpha**2)
    MidCont *= 12 * np.pi * (R**3) / (alpha**2)

    # Add results together
    result = InnerCont + OuterCont + MidCont
    return result

#############################################################################

'''
    Explicit chain rule conversions between angular derivatives
    in terms of spherical coordinates (theta, phi) and the 
    measure of angular separation (mu)
'''

# theta theta' derivative
def d2Y_dtdT(r, R, t, T, p, P, dLeg, ddLeg):
    dY = d2mu_dtdT(t, T, p, P) * dY_dmu(r, R, dLeg) + dmu_dt(t, T, p, P) * dmu_dT(t, T, p, P) * d2Y_ddmu(r, R, ddLeg)
    return dY
# theta phi' derivative
def d2Y_dtdP(r, R, t, T, p, P, dLeg, ddLeg):
    dY = d2mu_dtdP(t, T, p, P) * dY_dmu(r, R, dLeg) + dmu_dt(t, T, p, P) * dmu_dP(t, T, p, P) * d2Y_ddmu(r, R, ddLeg)
    return dY
# phi theta' derivative
def d2Y_dpdT(r, R, t, T, p, P, dLeg, ddLeg):
    dY = d2mu_dpdT(t, T, p, P) * dY_dmu(r, R, dLeg) + dmu_dp(t, T, p, P) * dmu_dT(t, T, p, P) * d2Y_ddmu(r, R, ddLeg)
    return dY
# phi phi' derivative
def d2Y_dpdP(r, R, t, T, p, P, dLeg, ddLeg):
    dY = d2mu_dpdP(t, T, p, P) * dY_dmu(r, R, dLeg) + dmu_dp(t, T, p, P) * dmu_dP(t, T, p, P) * d2Y_ddmu(r, R, ddLeg)
    return dY

#############################################################################

'''
    Linear combinations of derivatives of Upsilon make 
    up the covariance functions, as given in (B5)
'''

# Covariance of radial magnetic field measurements
#   C_{r r} (r, r')
def Crr(r, R, t, T, p, P):
    '''
        Inputs correspond to measurement positions:
            (r, theta, phi) = (r, t, p)
            (r', theta', phi') = (R, T, P)
    '''
    # Check for r <= R
    if (r > R):
        r, R = R, r
        t, T = T, t
        p, P = P, p
    # When theta = 0 exactly, theta and phi directions are
    #  undefined so this is avoided by setting theta to a small value
    #   [converges in all cases as theta -> 0]
    if (t == 0):
        t = 1e-8
    if (T == 0):
        T = 1e-8
    
    # Calculate mu
    mu = np.cos(t) * np.cos(T) + np.sin(t) * np.sin(T) * np.cos(p - P)
    # Negate computational error when incorrectly evaluating mu != +-1
    if (1 - mu < CutOff):
        mu = 1
    elif (1 + mu < CutOff):
        mu = -1
    
    # Calculate Legendre polynomial and derivatives of mu
    Pn, dPn, ddPn = LegP_finder(mu)
    
    # Finding the covariance (in nT^2)
    C = 1e18 * (mu0**2 / (r * R * 16 * R_P * np.pi**2)) * (2 * mu * dY_dmu(r, R, dPn) + (mu**2 - 1) * d2Y_ddmu(r, R, ddPn))
    # Output result
    return C

# Covariance of radial and colatitudinal magnetic field measurements
#   C_{r theta} (r, r')
def Crt(r, R, t, T, p, P):
    '''
        Inputs correspond to measurement positions:
            (r, theta, phi) = (r, t, p)
            (r', theta', phi') = (R, T, P)
    '''
    # Check for r <= R
    if (r > R):
        C = Ctr(R, r, T, t, P, p)
        return C
    else:
        # When theta = 0 exactly, theta and phi directions are
        #  undefined so this is avoided by setting theta to a small value
        #   [converges in all cases as theta -> 0]
        if (t == 0):
            t = 1e-8
        if (T == 0):
            T = 1e-8
            
        # Calculate mu
        mu = np.cos(t) * np.cos(T) + np.sin(t) * np.sin(T) * np.cos(p - P)
        # Negate computational error when incorrectly evaluating mu != +-1
        if (1 - mu < CutOff):
            mu = 1
        elif (1 + mu < CutOff):
            mu = -1
        
        # Calculate Legendre polynomial and derivatives of mu
        Pn, dPn, ddPn = LegP_finder(mu)

        # Finding the covariance (in nT^2)
        C = 1e18 * (mu0**2 / (r * R * 16 * R_P * np.pi**2)) * dmu_dT(t, T, p, P) * (dY_dmu(r, R, dPn) + R * d2Y_dRdmu(r, R, dPn))
        # Output result
        return C

# Covariance of radial and azimuthal magnetic field measurements
#   C_{r phi} (r, r')
def Crp(r, R, t, T, p, P):
    '''
        Same inputs and structure as Crt
    '''
    if (r > R):
        C = Cpr(R, r, T, t, P, p)
        return C
    else:
        if (t == 0):
            t = 1e-8
        if (T == 0):
            T = 1e-8
            
        mu = np.cos(t) * np.cos(T) + np.sin(t) * np.sin(T) * np.cos(p - P)
        if (1 - mu < CutOff):
            mu = 1
        elif (1 + mu < CutOff):
            mu = -1
        
        Pn, dPn, ddPn = LegP_finder(mu)
        
        C = 1e18 * (mu0**2 / (r * R * np.sin(T) * 16 * R_P * np.pi**2)) * dmu_dP(t, T, p, P) * (dY_dmu(r, R, dPn) + R * d2Y_dRdmu(r, R, dPn))
        return C

# Covariance of colatitudinal and radial magnetic field measurements
#   C_{theta r} (r, r')
def Ctr(r, R, t, T, p, P):
    '''
        Same inputs and structure as Crt
    '''
    if (r > R):
        C = Crt(R, r, T, t, P, p)
        return C
    else:
        if (t == 0):
            t = 1e-8
        if (T == 0):
            T = 1e-8
            
        mu = np.cos(t) * np.cos(T) + np.sin(t) * np.sin(T) * np.cos(p - P)
        if (1 - mu < CutOff):
            mu = 1
        elif (1 + mu < CutOff):
            mu = -1
        
        Pn, dPn, ddPn = LegP_finder(mu)
        
        C = 1e18 * (mu0**2 / (r * R * 16 * R_P * np.pi**2)) * dmu_dt(t, T, p, P) * (dY_dmu(r, R, dPn) + r * d2Y_drdmu(r, R, dPn))
        return C

#Covariance of colatitudinal magnetic field measurements
#   C_{theta theta} (r, r')
def Ctt(r, R, t, T, p, P):
    '''
        Same inputs and structure as Crr
    '''
    if (r > R):
        r, R = R, r
        t, T = T, t
        p, P = P, p
    if (t == 0):
        t = 1e-8
    if (T == 0):
        T = 1e-8

    mu = np.cos(t) * np.cos(T) + np.sin(t) * np.sin(T) * np.cos(p - P)
    if (1 - mu < CutOff):
        mu = 1
    elif (1 + mu < CutOff):
        mu = -1

    Pn, dPn, ddPn = LegP_finder(mu)

    C = 1e18 * (mu0**2 / (r * np.sin(t) * R * np.sin(T) * 16 * R_P * np.pi**2)) * (mu * d2Y_dpdP(r, R, t, T, p, P, dPn, ddPn) + r * dmu_dP(t, T, p, P) * dmu_dP(t, T, p, P) * d2Y_drdmu(r, R, dPn) + R * dmu_dp(t, T, p, P) * dmu_dp(t, T, p, P) * d2Y_dRdmu(r, R, dPn) + r * R * d2mu_dpdP(t, T, p, P) * d2Y_drdR(r, R, Pn))
    return C

#Covariance of colatitudinal and azimuthal magnetic field measurements
#   C_{theta phi} (r, r')
def Ctp(r, R, t, T, p, P):
    '''
        Same inputs and structure as Crt
    '''
    if (r > R):
        C = Cpt(R, r, T, t, P, p)
        return C
    else:
        if (t == 0):
            t = 1e-8
        if (T == 0):
            T = 1e-8
            
        mu = np.cos(t) * np.cos(T) + np.sin(t) * np.sin(T) * np.cos(p - P)
        if (1 - mu < CutOff):
            mu = 1
        elif (1 + mu < CutOff):
            mu = -1
        
        Pn, dPn, ddPn = LegP_finder(mu)
        
        C = 1e18 * (mu0**2 / (r * np.sin(t) * R * 16 * R_P * np.pi**2)) * (- mu * d2Y_dpdT(r, R, t, T, p, P, dPn, ddPn) + r * dmu_dp(t, T, p, P) * dmu_dT(t, T, p, P) * d2Y_drdmu(r, R, dPn) + R * dmu_dT(t, T, p, P) * dmu_dp(t, T, p, P) * d2Y_dRdmu(r, R, dPn) - r * R * d2mu_dpdT(t, T, p, P) * d2Y_drdR(r, R, Pn))
        return C

#Covariance of azimuthal and radial magnetic field measurements
#   C_{phi r} (r, r')
def Cpr(r, R, t, T, p, P):
    '''
        Same inputs and structure as Crt
    '''
    if (r > R):
        C = Crp(R, r, T, t, P, p)
        return C
    else:
        if (t == 0):
            t = 1e-8
        if (T == 0):
            T = 1e-8
            
        mu = np.cos(t) * np.cos(T) + np.sin(t) * np.sin(T) * np.cos(p - P)
        if (1 - mu < CutOff):
            mu = 1
        elif (1 + mu < CutOff):
            mu = -1
        
        Pn, dPn, ddPn = LegP_finder(mu)
        
        C = 1e18 * (mu0**2 / (r * np.sin(t) * R * 16 * R_P * np.pi**2)) * dmu_dp(t, T, p, P) * (dY_dmu(r, R, dPn) + r * d2Y_drdmu(r, R, dPn))
        return C

#Covariance of azimuthal and colatitudinal magnetic field measurements
#   C_{phi theta} (r, r')
def Cpt(r, R, t, T, p, P):
    '''
        Same inputs and structure as Crt
    '''
    if (r > R):
        C = Ctp(R, r, T, t, P, p)
        return C
    else:
        if (t == 0):
            t = 1e-8
        if (T == 0):
            T = 1e-8
            
        mu = np.cos(t) * np.cos(T) + np.sin(t) * np.sin(T) * np.cos(p - P)
        if (1 - mu < CutOff):
            mu = 1
        elif (1 + mu < CutOff):
            mu = -1

        Pn, dPn, ddPn = LegP_finder(mu)
        
        C = 1e18 * (mu0**2 / (r * R * np.sin(T) * 16 * R_P * np.pi**2)) * (- mu * d2Y_dtdP(r, R, t, T, p, P, dPn, ddPn) + R * dmu_dP(t, T, p, P) * dmu_dt(t, T, p, P) * d2Y_dRdmu(r, R, dPn) + r * dmu_dt(t, T, p, P) * dmu_dP(t, T, p, P) * d2Y_drdmu(r, R, dPn) - r * R * d2mu_dtdP(t, T, p, P) * d2Y_drdR(r, R, Pn))
        return C
    
#Covariance of azimuthal magnetic field measurements
#   C_{phi phi} (r, r')
def Cpp(r, R, t, T, p, P):
    '''
        Same inputs and structure as Crr
    '''
    if (r > R):
        r, R = R, r
        t, T = T, t
        p, P = P, p
    if (t == 0):
        t = 1e-8
    if (T == 0):
        T = 1e-8
    
    mu = np.cos(t) * np.cos(T) + np.sin(t) * np.sin(T) * np.cos(p - P)
    if (1 - mu < CutOff):
        mu = 1
    elif (1 + mu < CutOff):
        mu = -1
    
    Pn, dPn, ddPn = LegP_finder(mu)
    
    C = 1e18 * (mu0**2 / (r * R * 16 * R_P * np.pi**2)) * (mu * d2Y_dtdT(r, R, t, T, p, P, dPn, ddPn) - r * dmu_dt(t, T, p, P) * dmu_dT(t, T, p, P) * d2Y_drdmu(r, R, dPn) - R * dmu_dT(t, T, p, P) * dmu_dt(t, T, p, P) * d2Y_dRdmu(r, R, dPn) + r * R * d2mu_dtdT(t, T, p, P) * d2Y_drdR(r, R, Pn))
    return C
