#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A test program to learn python classes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import pylab
from numpy import *
from alm_classes import *

# Set up parameters

wavelength = 532e-9
n = 1.4607

R_ETM = 2312
R_primary = 4
R_secondary = -0.2

f_L3 = -0.2171
f_L4 = 1.0853
f_ETM = -R_ETM/(n-1)

L3toSecondary = 6.1549
SecondarytoPrimary = 1.9026
PrimarytoETM = 2.0

target_waist = 0.0115
target_waistpos = 2002.5 + L3toSecondary + SecondarytoPrimary + PrimarytoETM


# Initialize beam path object with wavelength, input beam parameters (defined by beam scan upstream from L3)

ALS_Path = beamPath()

ALS_Path.setWavelength(532e-9)

ALS_Path.setInputBeam(245e-6,-0.109)

# Calculate optimal position of in-air ALS L4
# Order of optics is L3 (z=0), L4, Secondary, Primary, ETM
# note current position of L4 is 0.944

L4_min_pos = 0.01
L4_max_pos = 2.0

L4_pos = arange(L4_min_pos,L4_max_pos,0.005)
waist_array = []
waistpos_array = []

for i in L4_pos:

    ALS_Path.clearElements()

    # Fixed elements: L3, secondary, primary, ETM

    Lens3 = beamComponent()
    Lens3.Lens(0.0,f_L3,'L3')
    ALS_Path.addElement(Lens3)

    Secondary = beamComponent()
    Secondary.Mirror(L3toSecondary,R_secondary,'Secondary')
    ALS_Path.addElement(Secondary)

    dist_sec2pri = beamComponent()
    dist_sec2pri.Dist(L3toSecondary+0.0001,SecondarytoPrimary,'Sec2Pri')
    ALS_Path.addElement(dist_sec2pri)

    Primary = beamComponent()
    Primary.Mirror(L3toSecondary+SecondarytoPrimary,R_primary,'Primary')
    ALS_Path.addElement(Primary)

    dist_pri2ETM = beamComponent()
    dist_pri2ETM.Dist(L3toSecondary+SecondarytoPrimary+0.0001,PrimarytoETM,'Pri2ETM')
    ALS_Path.addElement(dist_pri2ETM)

    ETM = beamComponent()
    ETM.Lens(L3toSecondary+SecondarytoPrimary+PrimarytoETM,f_ETM,'ETM')
    ALS_Path.addElement(ETM)


    # Variable elements: L4, distance from L3 to L4, distance from L4 to secondary

    dist_L3toL4 = beamComponent()
    dist_L3toL4.Dist(0.0001,i,'L3toL4')
    ALS_Path.addElement(dist_L3toL4)

    Lens4 = beamComponent()
    Lens4.Lens(i,f_L4,'L3')
    ALS_Path.addElement(Lens4)
    
    dist_L4toSec = beamComponent()
    dist_L4toSec.Dist(i+0.0001,L3toSecondary-i,'L4toSec')
    ALS_Path.addElement(dist_L4toSec)

    q_new = ALS_Path.propagateGaussianBeam()

    w,z = ALS_Path.unpackComplexBeamParameter(L3toSecondary+SecondarytoPrimary+PrimarytoETM,q_new)

    waist_array.append(w)
    waistpos_array.append(z)


from matplotlib.font_manager import FontProperties

fontP = FontProperties()
fontP.set_size('medium')

fignum=0
fignum=fignum+1
pylab.figure(fignum)

pylab.subplot(211)
pylab.plot(L4_pos,waist_array,'r--')
pylab.plot([0.0,2.0],[target_waist,target_waist],'k:')
pylab.plot(0.944,target_waist,'bs')

pylab.subplot(212)
pylab.plot(L4_pos,waistpos_array,'r--')
pylab.plot([0.0,2.0],[target_waistpos,target_waistpos],'k:')
pylab.plot(0.944,target_waistpos,'bs')


pylab.savefig('fit_ALS_telescope.png')


fignum=fignum+1
pylab.figure(fignum)
pylab.plot(waist_array,waistpos_array,'ko')
pylab.plot(target_waist,target_waistpos,'b*')

pylab.savefig('map_ALS_telescope.png')
