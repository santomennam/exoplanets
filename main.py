import pandas as pd
import math
import matplotlib
from matplotlib import pyplot as plt
from IPython.display import display

stefanBoltzmannsConstant = 5.67e-8
earthDistFromSun = 149.6e9  # km
earthEquilTemp = 278.5  # kelvin


class planetSystem:
    def __init__(self, name, luminosity, distance, stellarRadius, albedo, beta, radius, mass, emissivity=1):
        self.mass = mass #in earth masses
        self.radius = radius # earth radii
        self.name = name # name of planet
        self.luminosity = luminosity #in stellar units
        self.distance = distance #in au
        self.stellarRadius = stellarRadius #sun radii
        self.albedo = albedo
        self.beta = beta #fraction of the planet that re-radiates the absorbed flux, default 1
        self.emissivity = emissivity # usually ~1
        self.eTemp = equilibriumTemp(self.luminosity, self.stellarRadius, self.distance, self.albedo, self.beta)


def stellarFlux(luminosity, distance):
    if (type(luminosity == 'str')):
        luminosity = float(luminosity)
    return luminosity / (4 * math.pi * distance * distance)


# The factor beta is the fraction of the planet surface that re-radiates the absorbed flux, b = 1 for fast rotators and b ~= 0.5 for tidally locked planets without oceans or atmospheres
def equilibriumTemp(luminosity, distance, stellarRadius, albedo, beta, emissivity=1):
    return pow(((1 - albedo) * stellarFlux(luminosity, distance)) / (4 * beta * emissivity * stefanBoltzmannsConstant), 0.25)
    # see eqs 1 and 2 in Méndez, A. and Rivera-Valentín, E.G., 2017. The equilibrium temperature of planets in elliptical orbits. The Astrophysical Journal Letters, 837(1), p.L1.


df = pd.read_csv('exoData.csv')

#provided data is as follows:
# Planet	T_(eq)	T_(s)	A_(b)	G_(n)	kappa
# ------------------------------------------------
# Venus	     229	730	    0.760	0.990	2.213
# Earth	     255	288	    0.301	0.385	1.033
# Mars	     210	214	    0.250	0.073	0.948
# Titan	     84	94	0.265	0.338	1.027

#outputted T_s as of 6/22/22:
# earth: 287.558
# venus: 724.800
# mars: 213.7265
# titan: 92.7514
# (pretty good)

#two different formulas for the same thing, for sanity check.
def surfaceTemp(greenhouse, luminosity, beta, distance, albedo, emissivity=1):
    return earthEquilTemp * pow(((1 - albedo) * luminosity / (beta * emissivity *(1-greenhouse)* distance * distance)), 0.25)

def newSurfaceTemp(greenhouse, luminosity, beta, distance, albedo, emissivity=1):
    return earthEquilTemp * pow((1-albedo)/(1-greenhouse),0.25)*pow(stellarFlux(luminosity, distance)/stellarFlux(1,1),0.25)

def calc():
    df['calculatedEqEarthTemp'] = df.apply(lambda row: equilibriumTemp(row['luminosity'], row['planetSemimajor'], row['stellarRadius'], 0.301, 1), axis=1)
    df['calculatedEqVenusTemp'] = df.apply(lambda row: equilibriumTemp(row['luminosity'], row['planetSemimajor'], row['stellarRadius'], 0.760, 1), axis=1)
    df['calculatedEqMarsTemp'] = df.apply(lambda row: equilibriumTemp(row['luminosity'], row['planetSemimajor'], row['stellarRadius'], 0.250, 1), axis=1)
    df['calculatedEqTitanTemp'] = df.apply(lambda row: equilibriumTemp(row['luminosity'], row['planetSemimajor'], row['stellarRadius'], 0.265, 1), axis=1)

    df['calculatedSurfaceEarthTemp'] = df.apply(lambda row: surfaceTemp(0.385, row['luminosity'], 1, row['planetSemimajor'], 0.301, 1), axis=1)
    df['calculatedSurfaceVenusTemp'] = df.apply(lambda row: surfaceTemp(0.990, row['luminosity'], 1, row['planetSemimajor'], 0.760, 1), axis=1)
    df['calculatedSurfaceMarsTemp'] = df.apply(lambda row: surfaceTemp(0.073, row['luminosity'], 1, row['planetSemimajor'], 0.250, 1), axis=1)
    df['calculatedSurfaceTitanTemp'] = df.apply(lambda row: surfaceTemp(0.338, row['luminosity'], 1, row['planetSemimajor'], 0.265, 1), axis=1)


print("old earth: " + str(surfaceTemp(0.385,1,1,1,0.301,1))+"\nold venus: " + str(surfaceTemp(0.990,1,1,0.7233,0.760,1))+"\nold mars: " + str(surfaceTemp(0.073,1,1,1.5273,0.250,1))+"\nold titan: " + str(surfaceTemp(0.338,1,1,9.5,0.265,1)))
print("\nnew earth: " + str(newSurfaceTemp(0.385,1,1,1,0.301,1))+"\nnew venus: " + str(newSurfaceTemp(0.990,1,1,0.7233,0.760,1))+"\nnew mars: " + str(newSurfaceTemp(0.073,1,1,1.5273,0.250,1))+"\nnew titan: " + str(newSurfaceTemp(0.338,1,1,9.5,0.265,1)))

calc()
df.plot.scatter(x='calculatedSurfaceEarthTemp', y='calculatedEqEarthTemp', c='calculatedSurfaceEarthTemp', colormap='viridis')
plt.show()
#we get some pretty big values here. sus
