#
setwd("C:\\Users\\Xiang\\Google Drive\\Research\\RFIM_Simulation\\Simulation\\figs_code")

mfn = "..\\data\\MED_L64_T0.1_shi1.0_sHi0.32_r1000_sqT.csv"
med = read.csv(mfn, sep=',', skip=1)
with(med, plot(mc_step, mag) )
with(med, plot(mc_step, energy) )
with(med, plot(mc_step, domain_len) )

mfn = "..\\data\\MED_L128_T0.1_shi1.0_sHi0.32_r200_sqT.csv"
med = read.csv(mfn, sep=',', skip=1)
with(med, plot(mc_step, mag) )
with(med, plot(mc_step, energy) )
with(med, plot(mc_step, domain_len) )
