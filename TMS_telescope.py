#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# An optical beam path for the aLIGO Transmon Telescope IR path.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
from alm_classes import *

# set wavelength, index of refraction, lens focal lengths
#wavelength = 532e-9
#n = 1.4607

wavelength = 1064e-9
n = 1.449

R_ETM = 2245
R_primary = 4
R_secondary = -0.2

# Think that L101 is PLCX R=154.5mm, L102 is PLCC R=51.5mm

R_L101 = 0.1545
R_L102 = -0.0515

f_ETM = -R_ETM/(n-1)
f_L101 = R_L101/(n-1)
f_L102 = R_L102/(n-1)

TMS_Path = beamPath()

TMS_Path.setWavelength(wavelength)

# From LIGO-T0900385-v07, table 3.  ETMtoERM includes ETM thickness
ETMtoERM = 0.22
ERMtoPrimary = 1.65
Primary2Secondary = 1.9026
Secondary2L101 = 1.0974
L101toL102 = 0.246
L102toQPD1 = 0.439
QPD1toQPD2 = 0.369


# Calculate QPD sensitivity to beam displacements/angles
# Order of optics is ETM, ERM, primary, secondary, L101, L102, QPD1, QPD2

arm_waist = 0.012
arm_waistpos = 2160  # distance from ETM

dist_Waist2ETM = beamComponent()
dist_Waist2ETM.Dist(0.0,arm_waistpos,'Waist2ETM')
TMS_Path.addElement(dist_Waist2ETM)

ETM = beamComponent()
ETM.Lens(arm_waistpos,f_ETM,'ETM')
TMS_Path.addElement(ETM)

dist_ETM2ERM = beamComponent()
dist_ETM2ERM.Dist(arm_waistpos+0.0001,ETMtoERM,'ETM2ERM')
TMS_Path.addElement(dist_ETM2ERM)

ERM = beamComponent()
ERM.Slab(arm_waistpos+ETMtoERM,0.13,n,'ERM')
TMS_Path.addElement(ERM)

dist_ERM2Pri = beamComponent()
dist_ERM2Pri.Dist(arm_waistpos+ETMtoERM+0.0001,ERMtoPrimary,'ERMtoPrimary')
TMS_Path.addElement(dist_ERM2Pri)

Primary = beamComponent()
Primary.Mirror(arm_waistpos+ETMtoERM+ERMtoPrimary,R_primary,'Primary')
TMS_Path.addElement(Primary)

dist_pri2sec = beamComponent()
dist_pri2sec.Dist(arm_waistpos+ETMtoERM+ERMtoPrimary+0.0001,Primary2Secondary,'Pri2Sec')
TMS_Path.addElement(dist_pri2sec)

Secondary = beamComponent()
Secondary.Mirror(arm_waistpos+ETMtoERM+ERMtoPrimary+Primary2Secondary,R_secondary,'Secondary')
TMS_Path.addElement(Secondary)

dist_sec2L101 = beamComponent()
dist_sec2L101.Dist(arm_waistpos+ETMtoERM+ERMtoPrimary+Primary2Secondary+0.0001,Secondary2L101,'Sec2L101')
TMS_Path.addElement(dist_sec2L101)

L101 = beamComponent()
L101.Lens(arm_waistpos+ETMtoERM+ERMtoPrimary+Primary2Secondary+Secondary2L101,f_L101,'L101')
TMS_Path.addElement(L101)

dist_L101toL102 = beamComponent()
dist_L101toL102.Dist(arm_waistpos+ETMtoERM+ERMtoPrimary+Primary2Secondary+Secondary2L101+0.0001,L101toL102,'L101toL102')
TMS_Path.addElement(dist_L101toL102)

L102 = beamComponent()
L102.Lens(arm_waistpos+ETMtoERM+ERMtoPrimary+Primary2Secondary+Secondary2L101+L101toL102,f_L102,'L102')
TMS_Path.addElement(L102)

dist_L102toQPD1 = beamComponent()
dist_L102toQPD1.Dist(arm_waistpos+ETMtoERM+ERMtoPrimary+Primary2Secondary+Secondary2L101+L101toL102+0.0001,L102toQPD1,'L102toQPD1')
TMS_Path.addElement(dist_L102toQPD1)
"""
dist_QPD1toQPD2 = beamComponent()
dist_QPD1toQPD2.Dist(arm_waistpos+ETMtoERM+ERMtoPrimary+Primary2Secondary+Secondary2L101+L101toL102+L102toQPD1,QPD1toQPD2,'QPD1toQPD2')
TMS_Path.addElement(dist_QPD1toQPD2)

"""
TMS_Path.setInputBeam(arm_waist,0.0)

q_new = TMS_Path.propagateGaussianBeam()

#TMS_Path.printElementLabels()

print 'Complex beam parameter is ', q_new
print

w,z = TMS_Path.unpackComplexBeamParameter(arm_waistpos+ETMtoERM+ERMtoPrimary+Primary2Secondary+Secondary2L101+L101toL102+L102toQPD1,q_new)
#w,z = TMS_Path.unpackComplexBeamParameter(arm_waistpos+ETMtoERM+ERMtoPrimary+Primary2Secondary+Secondary2L101+L101toL102+L102toQPD1+QPD1toQPD2,q_new)
#w,z = TMS_Path.unpackComplexBeamParameter(arm_waistpos+ETMtoERM+ERMtoPrimary+Primary2Secondary+Secondary2L101,q_new)

#print 'Current beam position is ', ETMtoERM+ERMtoPrimary+Primary2Secondary+Secondary2L101+L101toL102+L102toQPD1

# Returns the 1/e^2 radius (in power) of the beam at position z1 for the complex beam parameter q defined at z2
z1 = arm_waistpos+ETMtoERM+ERMtoPrimary+Primary2Secondary+Secondary2L101+L101toL102+L102toQPD1+QPD1toQPD2
z2 = arm_waistpos+ETMtoERM+ERMtoPrimary+Primary2Secondary+Secondary2L101+L101toL102+L102toQPD1

#print z1-arm_waistpos, z2-arm_waistpos, TMS_Path.getLastElementPos()-arm_waistpos

print 'QPD1 position is ', ETMtoERM+ERMtoPrimary+Primary2Secondary+Secondary2L101+L101toL102+L102toQPD1
print 'QPD2 position is ', ETMtoERM+ERMtoPrimary+Primary2Secondary+Secondary2L101+L101toL102+L102toQPD1+QPD1toQPD2
print

z1 = TMS_Path.getLastElementPos()
beam_size = TMS_Path.getBeamSize(z1,q_new,z1)

print 'Waist size, location: ', w, z-arm_waistpos
print 'Position of last element is ', z1-arm_waistpos
print 'Beam size is ', beam_size
print

print 'Angle sensitivity is ', TMS_Path.getAngleSensitivity()/w
print 'Lateral sensitivity is ', TMS_Path.getLateralSensitivity()/w
print
