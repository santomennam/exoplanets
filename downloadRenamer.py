import sys
import pandas as pd

columnDict = {"pl_name":"Name","hostname":"starName","default_flag":"defaultPlanetEntry","sy_snum":"numStars","sy_pnum":"numPlanets","soltype":"solution","pl_controv_flag":"controversial",
             "pl_orbper":"orbitalPeriod","pl_orbsmax":"planetSemimajor","pl_rade":"planetEarthRads","pl_masse":"earthMass",
              "pl_bmasse":"bestEarthMass","pl_bmassprov":"bestMassProvenance","pl_orbeccen":"eccentricity","pl_eqt":"equilTemp",
              "st_teff":"stellarEffTemp","st_rad":"stellarRadius","st_met":"stellarMetallicity","st_metratio":"stellarMetallicityRatio","st_lum":"luminosity","st_age":"stellarAge","sy_dist":"distance" }

path = sys.argv[1]
saveName = sys.argv[2]
df = pd.read_csv(path)
print("Columns before:")
print(df.columns)
df.rename(columns=columnDict,inplace=True)
print("Columns after rename:")
print(df.columns,"\n")
print("df:")
print(df)

df.to_csv(saveName)

#pl_name:        Planet Name
#hostname:       Host Name
#default_flag:   Default Parameter Set
#sy_snum:        Number of Stars
#sy_pnum:        Number of Planets
#soltype:        Solution Type
#pl_controv_flag: Controversial Flag
#pl_orbper:      Orbital Period [days]
#pl_orbsmax:     Orbit Semi-Major Axis [au])
#pl_rade:        Planet Radius [Earth Radius]
#pl_masse:       Planet Mass [Earth Mass]
#pl_bmasse:      Planet Mass or Mass*sin(i) [Earth Mass]
#pl_bmassprov:   Planet Mass or Mass*sin(i) Provenance
#pl_orbeccen:    Eccentricity
#pl_eqt:         Equilibrium Temperature [K]
#st_teff:        Stellar Effective Temperature [K]
#st_rad:         Stellar Radius [Solar Radius]
#st_met:         Stellar Metallicity [dex]
#st_metratio:    Stellar Metallicity Ratio
#st_lum:         Stellar Luminosity [log(Solar)]
#st_age:         Stellar Age [Gyr]
#sy_dist:        Distance [pc]
