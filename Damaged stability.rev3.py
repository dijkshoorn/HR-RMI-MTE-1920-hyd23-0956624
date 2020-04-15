# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 12:51:18 2020

@author: mauri
"""

"PYTHON LIBRARIES"
import numpy as np

"VESSEL A834"
T = 5.5 #[m], Draft (Mean) - From PIAS;Hydrotables (Extensive)
Nabla = 80 * 18 * 5.5 #[], Deplacemnet - From PIAS;Hydrotables (Extensive)
A_wl = 80 * 18 #[m^2], Waterplane area - From PIAS;Hydrotables (Extensive)
CoF = 40 #[m], Center of flotation - From ord. 0

MoI_t = (1 / 12) * 80 * 18**3 #[m^4], Transverse moment of inertia of waterplane - From PIAS;Hydrotables (Extensive)
MoI_L = (1 / 12) * 18 * 80**3 #[m^4], Longitudinal moment of inertia of waterplane - From PIAS;Hydrotables (Extensive)

MoI_t_c1 = (1 / 12) * 16 * 9**3 #[m^4], Transverse moment of inertia ()
MoI_L_c1 = (1 / 12) * 9 * 16**3 #[m^4], Longitudinal moment of inertia ()


"PERMEABILITY"
Mu_v = 0.8 #[-], Volume permeability  
Mu_a = 1.0 #[-], Area permeability

"COMPARTMENTS"
V_c1 = 16 * 9 * 5.5 #[m^3], Volume damaged compartment

VCG_c1 = 2.75 #[m], Vertical center of gravity from base
TCG_c1 = 4.5 #[m], Transverse center of gravity 
LCG_c1 = (80 / 2) - (16 / 2) + (80 / 2) #[m], Longitudinal center of gravity 

V_ct = V_c1 
A_c = 16 * 9 #[m^2], Damaged compartment area - From PIAS;Tanktable Areas 

"STABILITY"
KB_0 = 2.75 #[m], Vertical centre of bouyancy (undamaged) - From PIAS;Hydrotables (Extensive) 
KG = 4.5 #[m], Vertical centre of gravity (undamaged) - From PIAS;Loading;Weightlist

"CALCULATIONS"
"NEW T"
Delta_T = (Mu_v * V_ct) / (A_wl - (Mu_v * A_c) ) #[m], Delta draft

"VERTICAL MOVEMENT OF CENTER OF BUOYANCY"
VCG_ct = VCG_c1
h = T + (0.5 * Delta_T) - VCG_ct #[m], New vertical center of gravity
ZB = Mu_v * V_ct * h / Nabla

"TRANSVERSE MOVEMENT OF CENTER OF BUOYANCY"
c = (V_c1 * TCG_c1) / V_ct #[m], Transverse center of gravity of combined compartments
y_1 = (Mu_a * A_c * c) / (A_wl - Mu_a * A_c)
b = y_1 + c
y_B = Mu_v * V_ct * b / Nabla

"LONGITUDINAL MOVEMENT OF CENTER OF BUOYANCY"
a = ( (V_c1 * LCG_c1) / V_ct) - CoF #[m], Distance between center of flotation & Longitudinal center of gravity of combined compartments
x_1 = (Mu_a * A_c * a) / (A_wl - Mu_a * A_c)
l = x_1 + a
x_B = (Mu_v * V_ct * l) / Nabla

"CALCULATION OF NEW STABILITY PARAMETERS"
"KB_1"
KB_1 = KB_0 + ZB

"BM VALUES"
MoI_t_ct = MoI_t_c1
MoI_L_ct = MoI_L_c1

I_t1 = MoI_t - (Mu_a * MoI_t_ct) - (Mu_a * A_c * c**2) - (A_wl - Mu_a * A_c) * y_1**2
I_L1 = MoI_L - (Mu_a * MoI_L_ct) - (Mu_a * A_c * a**2) - (A_wl - Mu_a * A_c) * x_1**2

"BM_1"
BM_1 = I_t1 / Nabla

"BM_L1"
BM_L1 = I_L1 / Nabla

"GM VALUES"
"GM_1"
GM_1 = KB_1 + BM_1 - KG

"GM_L1"
GM_L1 = KB_1 + BM_L1 - KG

"ANGLES"
tan_Phi = y_B / GM_1
tan_Theta = x_B / GM_L1

Phi = np.degrees(np.arctan(tan_Phi)) #[degrees], Capsizing angle
Theta = np.degrees(np.arctan(tan_Theta)) #[degrees], Trim angle

"RESULTS"
print("Delta T =", Delta_T,
      "h =", h,
      "ZB =", ZB,
      "c =", c,
      "y1 =", y_1,
      "b =", b,
      "yB =", y_B,
      "a =", a,
      "x1 =", x_1,
      "l =", l,
      "xB =", x_B,
      "KB1 =", KB_1,
      "It1 =", I_t1,
      "Il1 =", I_L1,
      "BM1 =", BM_1,
      "BML1 =", BM_L1,
      "GM1 =", GM_1,
      "GML1 =", GM_L1,
      "Phi =", Phi,
      "Theta =", Theta)