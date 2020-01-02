import pysm
from pysm.nominal import models
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib as mpl
from pysm.common import convert_units
from NPTFit import create_mask as cm
import matplotlib.style
mpl.style.use('classic')

nside = 128
nu = np.logspace(np.log10(30), np.log10(400), 30); Nf = len(nu)
coefficients = convert_units("uK_RJ", "uK_CMB", nu)

def convert_unit(map):
    for i in range(0,Nf):
        map[i] = map[i]*coefficients[i]
    return map

def B_map(maps):
    alms = hp.map2alm(maps)
    Bmap = hp.alm2map(alms[2], nside = nside, verbose = False)
    return Bmap

s_rms = []; d_rms = []; t_rms = [];Cl = []; a_rms = []; c_rms = []
for i in range(1):
    for s in range(1,2): #3
        s1_config = models("s%s"%s, nside)

        for d in range(1,3): #7
            d1_config = models("d%s"%d, nside)

            for a in range(2,3): #2
                a1_config = models("a%s"%a, nside)
                c1_config = models("c1", 128)

                sky_config = {
                                'synchrotron' : s1_config,
                                'dust' : d1_config,
#                                 'freefree' : f1_config,
                                'cmb' : models("c1", nside),
                                'ame' : a1_config,
                            }
                sky = pysm.Sky(sky_config);
                synchrotron = (convert_unit(sky.synchrotron(nu)));
                dust = (convert_unit(sky.dust(nu)))
                ame = convert_unit(sky.ame(nu))
                cmb = convert_unit(sky.cmb(nu))
                total = (convert_unit(sky.signal()(nu)) - cmb)
    
                total_Bmap = []; dust_Bmap = []; sync_Bmap = []; ame_Bmap = []; cmb_Bmap = []  
                for i in range (len(nu)):
                    total_Bmap.append(B_map(total[i]))
                    dust_Bmap.append(B_map(dust[i]))
                    sync_Bmap.append(B_map(synchrotron[i]))
                    ame_Bmap.append(B_map(ame[i]))
                    cmb_Bmap.append(B_map(cmb[i]))
                    
                
                total_fgnd = np.std(total_Bmap, axis = 1); t_rms.append(total_fgnd)
                dust_rms = np.std(dust_Bmap,axis = 1); d_rms.append(dust_rms)
                sync_rms = np.std(sync_Bmap,axis = 1); s_rms.append(sync_rms)
                ame_rms = np.std(ame_Bmap, axis = 1);  a_rms.append(ame_rms)
                cmb_rms = np.std(cmb_Bmap, axis = 1); c_rms.append(ame_rms)
np.savetxt('./dust_rms.txt', d_rms)
np.savetxt('./sync_rms.txt', s_rms)
np.savetxt('./ame_rms.txt', a_rms)
np.savetxt('./cmb_rms.txt', c_rms)
np.savetxt('./total_rms.txt', t_rms)
                
#plt.figure(figsize = (10, 8))
#plt.loglog(nu, total_fgnd[:], label = "Total fgnd", lw=3)
#plt.loglog(nu, dust_rms[:], label = "dust", lw=2)
#plt.loglog(nu, sync_rms[:], label = "synchrotron", lw=2)
#plt.loglog(nu, ame_rms[:],'purple' ,label = "ame", lw=2)
#plt.loglog(nu, cmb_rms[:],'k' ,label = "CMB", lw=2)


#plt.xlabel(r"$\nu$ (GHz)", fontsize = 20)
#plt.ylabel(r"$B_{\rm RMS}$ ($\mu K$)", fontsize = 20)
#plt.legend(loc = 'upper center',ncol=3)
#plt.ylim([0.001, 0.5*1000]); plt.xlim(30,410)
#ax=plt.gca();
#ax.spines['bottom'].set_linewidth(2);
#ax.spines['left'].set_linewidth(2);
#ax.spines['right'].set_linewidth(2);
#ax.spines['top'].set_linewidth(2);
#plt.show()