#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# An optical model of the aLIGO TMS telescope to calculate ray tracing
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
from scipy import constants
from alm_classes import *

# set wavelength, index of refraction, lens focal lengths
wavelength = 532e-9
n = 1.4607

#wavelength = 1064e-9
#n = 1.449

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

# From LIGO-T0900385-v07, table 3.
ETMtoERM = 0.02
ERM_thickness = 0.13
ERMtoPrimary = 1.65
Primary2Secondary = 1.9026
Secondary2L101 = 1.0974
L101toL102 = 0.246
L102toQPD1 = 0.439
QPD1toQPD2 = 0.369


# Assume beams start at normal incidence to ETM HR surface, and acquire an angle relative to the 
# normal axis when leaving the ETM AR surface (which is vertically wedged at 0.07deg)

ETM_wedge = 0.07 * constants.pi / 180  # convert wedge angle to radians

# angle of beam with respect to normal is sin(ETM_wedge) / sin(theta) = n1/n2, 
# where n1 = 1 for vacuum and n2 = index of ref. for fused silica at whatever wavelength
# assume small angle approximation holds, sin(theta)=theta

# the angle is so small that travel from ETM AR surface to ERM back surface does not add any appreciable displacement (z~0.35m),
# so starting z is back surface of ERM

start_angle = n*sin(ETM_wedge)
start_pos = 0.0

dist_ERM2Pri = beamComponent()
dist_ERM2Pri.Dist(0.0,ERMtoPrimary,'ERMtoPrimary')
TMS_Path.addElement(dist_ERM2Pri)

Primary = beamComponent()
Primary.Mirror(ERMtoPrimary,R_primary,'Primary')
TMS_Path.addElement(Primary)

dist_pri2sec = beamComponent()
dist_pri2sec.Dist(ERMtoPrimary+0.0001,Primary2Secondary,'Pri2Sec')
TMS_Path.addElement(dist_pri2sec)

Secondary = beamComponent()
Secondary.Mirror(ERMtoPrimary+Primary2Secondary,R_secondary,'Secondary')
TMS_Path.addElement(Secondary)

dist_sec2L101 = beamComponent()
dist_sec2L101.Dist(ERMtoPrimary+Primary2Secondary+0.0001,Secondary2L101,'Sec2L101')
TMS_Path.addElement(dist_sec2L101)

L101 = beamComponent()
L101.Lens(ERMtoPrimary+Primary2Secondary+Secondary2L101,f_L101,'L101')
TMS_Path.addElement(L101)

dist_L101toL102 = beamComponent()
dist_L101toL102.Dist(ERMtoPrimary+Primary2Secondary+Secondary2L101+0.0001,L101toL102,'L101toL102')
TMS_Path.addElement(dist_L101toL102)

L102 = beamComponent()
L102.Lens(ERMtoPrimary+Primary2Secondary+Secondary2L101+L101toL102,f_L102,'L102')
TMS_Path.addElement(L102)

dist_L102toQPD1 = beamComponent()
dist_L102toQPD1.Dist(ERMtoPrimary+Primary2Secondary+Secondary2L101+L101toL102+0.0001,L102toQPD1,'L102toQPD1')
TMS_Path.addElement(dist_L102toQPD1)

dist_QPD1toQPD2 = beamComponent()
dist_QPD1toQPD2.Dist(ERMtoPrimary+Primary2Secondary+Secondary2L101+L101toL102+L102toQPD1,QPD1toQPD2,'QPD1toQPD2')
TMS_Path.addElement(dist_QPD1toQPD2)
"""
"""
TMS_Path.setInputBeamPointing(start_pos,start_angle)

#TMS_Path.printElementLabels()

TMS_Path.traceRay()

#print 'QPD1 position is ', ERMtoPrimary+Primary2Secondary+Secondary2L101+L101toL102+L102toQPD1
#print 'QPD2 position is ', ERMtoPrimary+Primary2Secondary+Secondary2L101+L101toL102+L102toQPD1+QPD1toQPD2
#print

z1 = TMS_Path.getLastElementPos()

print 'Last element is ', TMS_Path.getLastElementLabel()

print 'Ray displacement at last element is ', TMS_Path.final_x
print 'Angle at last element is ', TMS_Path.final_theta
print


