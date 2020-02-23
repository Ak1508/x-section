import numpy as np
import ROOT as R
from array import array
import time, string

startTime = time.time()
# =====================
#     some constants 
# =====================

mp     = 0.9382723
mp2    = mp*mp
radcon = 180./3.141592654
sigave = 0.0
#============      ====================
#  mc_input.dat file 
#============      ====================

maxev   =  500000
dxp     =  60
dyp     =  65
delup   =  30
deldown = -20
cs_falg =  0
rc_flag =  1

# =======================================
#          read run info 
# =======================================
print ("We need run info" )

targetName = input ("Enter the target name: ")

tar_name_id ={'carbon': 4,
              'ld2'   : 15,
              'lh2'   : 11}

target = tar_name_id[targetName]

#print (target)

#print (type(target))


#target     = int(input(" Which target id you want to put? "))
#beamEnergy = float(input(" What is your beam energy 10.602? "))



beamEnergy = 10.602
pcen       = float(input(" What is your Spectrometer Momentum in GeV/C? "))
hsec       = round(pcen - 0.017 * pcen, 3) # spectrometer off by 1.7% source Abel ???
#thetac     = float(input(" What is your Spectrometer angle? " ))
thetac     = 21.0
#print (hsec)
bmcur1 = bmcur2 = bcm1charge = bcm2charge = cltime = eltime = trackeff = trigeff = rate = 1
bcmavecharge = (bcm1charge+bcm2charge)/2.         # !!! Average over BCMs !!!
bmcur        = (bmcur1+bmcur2)/2.

thetacrad    = thetac/radcon # converting central angle to radian 


#print (thetacrad)
dep          = (delup-deldown)/100.*hsec


#==================  ===============
# Read in target data from file 
#=================  ================
tar_id         = 0
density        = 0 
tar_length     = 0
at_mass        = 0
at_no          = 0
aerial_density = 0

with open ('/w/hallc-scifs17exp/xem2/abishek/monte-carlo/mc-reweight/input/target/targetdata_2017.dat', 'r') as f:
    for line in f:
        data = line.split()
        #print (type(data[0]))
        if int(data[0]) == int(target):

            tar_id         = int(data[0])
            density        = float(data[1])
            tar_length     = float(data[2])
            at_mass        = float(data[3])
            at_no          = float(data[4])
            aerial_density = float(data[5])

            #print (data[0], data[1], data[2], data[3], data[4], data[5])
outputName = input(" What is your MC RootFile Name? ")
#======================================              ===============================================
#     laoding output file from rc-extrnal to calculate radiative correction by linear interpolation #                                 
#======================================              ================================================
ebeam, wsqr, theta, x_sec, radCorrfactor = np.loadtxt('/w/hallc-scifs17exp/xem2/abishek/monte-carlo/rc-externals/output/rad-corr-data/%s_21_deg.dat' %(outputName), delimiter = '\t', unpack = True)


xt = R.TGraph2D() # here I creating Object for x-sec 
rt = R.TGraph2D() # here I creating Object for radiactive-correc

# def interpolation():
# I'm looping over all the rows in table and assigning to the object created above
for i in range(len(ebeam)):
    xt.SetPoint(i, wsqr[i], theta[i], x_sec[i])
    rt.SetPoint(i, wsqr[i], theta[i], radCorrfactor[i])

lumdata = (aerial_density*6.022137e-10/at_mass)*(bcmavecharge/1.602177e-13)

# =============================================
#        Reading the Monte-Carlo rootfile    #
# =============================================

fName = "/w/hallc-scifs17exp/xem2/abishek/monte-carlo/mc-single-arm/worksim/%s_%s.root" %(outputName,str(pcen).replace('.','p'))

print ("your output file name is: %s_%s.root " %(outputName,str(pcen).replace('.','p')))

f        = R.TFile(fName,"READ")
t        = f.Get("h1411")
nentries = t.GetEntries()

print (nentries)

#=========================================
#     output Root Tree                 #
#================ ========================

file = R.TFile("%s_%s.root" %(outputName,str(pcen).replace('.','p')),"RECREATE") ## this is output rootfile
tree = R.TTree("tree", "tree")

xfoc       = array('d', [0])
yfoc       = array('d', [0])
dxdz       = array('d', [0])
dydz       = array('d', [0])
ztarini    = array('d', [0])
ytarini    = array('d', [0])
delini     = array('d', [0])
yptarini   = array('d', [0])
xptarini   = array('d', [0])
ztar       = array('d', [0])
ytarrec    = array('d', [0])
delrec     = array('d', [0])
yptarrec   = array('d', [0])
xptarrec   = array('d', [0])
xtarini    = array('d', [0])
fail_id    = array('d', [0])
born       = array('d', [0])
rci        = array('d', [0])
#rce        = array('d', [0])
#csback     = array('d', [0])
hse        = array('d', [0])
hstheta    = array('d', [0])
#sigmac     = array('d', [0])
q2         = array('d', [0])
w2         = array('d', [0])


#=======================    ================
#   Setting branch for output rootFile
#======================    =================


tree.Branch('xfoc',     xfoc,     'xfoc/D')
tree.Branch('yfoc',     yfoc,     'yfoc/D')
tree.Branch('dxdz',     dxdz,     'dxdz/D')
tree.Branch('dydz',     dydz,     'dydz/D')
tree.Branch('ztarini',  ztarini,  'ztarini/D')
tree.Branch('yini',     ytarini,  'ytarini/D')
tree.Branch('dppi',     delini,   'delini/D')
tree.Branch('yptarini', yptarini, 'yptarini/D')
tree.Branch('xptarini', xptarini, 'xptarini/D')
tree.Branch('zrec',     ztar,     'ztar/D')
tree.Branch('yrec',     ytarrec,  'ytarrec/D')
tree.Branch('dppr',     delrec,   'delrec/D')
tree.Branch('yprec',    yptarrec, 'yptarrec/D')
tree.Branch('xprec',    xptarrec, 'xptarrec/D')
tree.Branch('xtarini',  xtarini,  'xtarini/D')
tree.Branch('fail_id',  fail_id,  'fail_id/D')
tree.Branch('born',     born,     'born/D')
tree.Branch('rci',      rci,      'rci/D')
#tree.Branch('rce',   rce, 'rce/D')
#tree.Branch('csback',csback, 'csback/D')
tree.Branch('hse',      hse,      'hse/D')
tree.Branch('hstheta',  hstheta,  'hstheta/D')
#tree.Branch('sigmac', sigmac, 'sigmac/D')
tree.Branch('q2',       q2,       'q2/D')
tree.Branch('w2',       w2,       'w2/D')
tree.Branch('born',     born,     'born/D')


sigav = 0
ngen  = 0

for entry in range(nentries):
    if entry       ==10000: print (' analyzed 10000 events')
    if entry       ==50000: print (' analyzed 50000 events')
    if entry % 100000 == 0: print (' analyzed', entry,'events')

    t.GetEntry(entry)

    lpsxfp     = t.GetLeaf('psxfp').GetValue(0)
    lpsyfp     = t.GetLeaf('psyfp').GetValue(0)
    lpsxpfp    = t.GetLeaf('psxpfp').GetValue(0)
    lpsypfp    = t.GetLeaf('psypfp').GetValue(0)
    lpsztari   = t.GetLeaf('psztari').GetValue(0)
    lpsytari   = t.GetLeaf('psytari').GetValue(0)
    lpsdeltai  = t.GetLeaf('psdeltai').GetValue(0)
    lpsyptari  = t.GetLeaf('psyptari').GetValue(0)
    lpsxptari  = t.GetLeaf('psxptari').GetValue(0)
    lpsztar    = t.GetLeaf('psztar').GetValue(0)
    lpsytar    = t.GetLeaf('psytar').GetValue(0)
    lpsdelta   = t.GetLeaf('psdelta').GetValue(0)
    lpsxptar   = t.GetLeaf('psxptar').GetValue(0)
    lpsyptar   = t.GetLeaf('psyptar').GetValue(0)
    lpsxtari   = t.GetLeaf('psxtari').GetValue(0)
    lfry       = t.GetLeaf('fry').GetValue(0)
    lxsnum     = t.GetLeaf('xsnum').GetValue(0)
    lysnum     = t.GetLeaf('ysnum').GetValue(0)
    lxsieve    = t.GetLeaf('xsieve').GetValue(0)
    lysieve    = t.GetLeaf('ysieve').GetValue(0)
    lstop_id   = t.GetLeaf('stop_id').GetValue(0)

    xfoc[0]    = lpsxfp
    yfoc[0]    = lpsyfp
    dxdz[0]    = lpsxpfp
    dydz[0]    = lpsypfp
    ztarini[0] = lpsztari
    ytarini[0] = lpsytari
    delini[0]  = lpsdeltai
    yptarini[0]= lpsyptari
    xptarini[0]= lpsxptari
    ztar[0]    = lpsztar
    ytarrec[0] = lpsytar
    delrec[0]  = lpsdelta
    yptarrec[0]= lpsyptar
    xptarrec[0]= lpsxptar
    xtarini[0] = lpsxtari
    fail_id[0] = lstop_id


    hse[0]      = hsec*(1.+ delini[0]/100.)

    ypcor       = 0.0
    
    yptarrec[0] = yptarrec[0]-ypcor 
    
    thetaini    = np.arccos(np.cos(thetacrad+yptarini[0])*np.cos(xptarini[0])) 

    hstheta[0]  = np.arccos(np.cos(thetacrad+yptarrec[0])*np.cos(xptarrec[0]))

    sin2        = np.power(np.sin(thetaini/2.), 2) 

    nu          = beamEnergy - hse[0]
       
    q2[0]       = 4.0*hse[0]*beamEnergy*sin2

    w2[0]       = mp2 + 2.*mp*nu-q2[0]

    eff_cer     = 1 # 0.998 Replacing by 1 coz we r correcting our data

    eff_cal     = 1 #0.96503+0.75590e-01*hse[0]-0.65283e-01*hse[0]**2 + 0.26938e-01*hse[0]**3-0.53013e-02*hse[0]**4+0.39896e-03*hse[0]**5

    dt          = thetaini - thetacrad

    phasespcor  = 1./np.power(np.cos(dt), 3)
    
    born_val    = xt.Interpolate(w2[0] ,thetaini*radcon)

    born[0]     = xt.Interpolate(w2[0] ,thetaini*radcon)*phasespcor*eff_cer*eff_cal
     
    rci[0]      = rt.Interpolate(w2[0] ,thetaini*radcon)

    tree.Fill()

    if abs(xptarini[0])< dxp and abs(yptarini[0])< dyp and delini[0] > deldown and delini[0] < delup and born_val >= 0:
        sigave  = sigave + born_val  
        ngen    = ngen + 1

file.Write()
file.Close()

sigave          = sigave/ngen
phase_space     = 4.0*dxp*dyp*dep/1000.0
sigtot          = sigave*phase_space 
lummc           = ngen/sigtot*1000.0
lumfract        = lumdata/lummc
fract           = lumdata*phase_space/ngen/1000.00 

print (" Scale factor for MC ", fract)

with open("%s_%s.txt" %(outputName, str(pcen).replace('.','p')), 'w') as fout:
    fout.write(str(fract))
          

print ('\nThe analysis took %.3f minutes\n' % ((time.time() - startTime) / (60.)))  

