import pandas as pd
import math
import matplotlib
from matplotlib import pyplot as plt
from IPython.display import display

stefanBoltzmannsConstant = 5.67e-8
earthDistFromSun = 149.6e9  # km
earthEquilTemp = 278.5  # kelvin

def stellarFlux(luminosity, distance):
    if (type(luminosity == 'str')):
        luminosity = float(luminosity)
    return luminosity / (4 * math.pi * distance * distance)


# The factor beta is the fraction of the planet surface that re-radiates the absorbed flux, b = 1 for fast rotators and b ~= 0.5 for tidally locked planets without oceans or atmospheres
def equilibriumTemp(luminosity, distance, stellarRadius, albedo, beta, emissivity=1):
    return pow(((1 - albedo) * stellarFlux(luminosity, distance)) / (4 * beta * emissivity * stefanBoltzmannsConstant), 0.25)
    # see eqs 1 and 2 in Méndez, A. and Rivera-Valentín, E.G., 2017. The equilibrium temperature of planets in elliptical orbits. The Astrophysical Journal Letters, 837(1), p.L1.

def surfaceTemp(greenhouse, luminosity, beta, distance, albedo, emissivity=1):
    return earthEquilTemp * pow(((1 - albedo) * luminosity / (beta * emissivity * (1 - greenhouse) * distance * distance)), 0.25)


def newSurfaceTemp(greenhouse, luminosity, beta, distance, albedo, emissivity=1):
    return earthEquilTemp * pow((1 - albedo) / (1 - greenhouse), 0.25) * pow(stellarFlux(luminosity, distance) / stellarFlux(1, 1), 0.25)


class similarityProfile:
    def __init__(self, weightExponents):
        self.weightExponents = weightExponents  # list of list of exponents
        # list of exponents for: mass, radius, eTemp, calculatedEqEarthTemp, calculatedEqVenusTemp, calculatedEqMarsTemp, calculatedEqTitanTemp, calculatedSurfaceEarthTemp,
        # calculatedSurfaceVenusTemp, calculatedSurfaceMarsTemp, calculatedSurfaceTitanTemp

        # sample exponents:
        # earth:         [1,1,1,1,0,0,0,1,0,0,0]
        # venus:         [1,1,1,0,1,0,0,0,1,0,0]
        # mars:          [1,1,1,0,0,1,0,0,0,1,0]
        # titan:         [1,1,1,0,0,0,1,0,0,0,1]
        # no temp:       [1,1,1,0,0,0,0,0,0,0,0]
        # body-neutral:  [1,1,1,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]


#  SI = product(i=1->m) of (1-abs((xp_i-x_i)/(xp_i+x_i)))^w_i
# where xp is the planetary reference value, x_i is the actual value of the comparison planet, and m is the number of properties.
# takes a list of planets, each with a list of list of exponents and their own characteristics. these are the base planets.
# also takes a dataframe row and compares the reference planets to the real database planet.
def similarityIndex(referencePlanetSystems, row):
    if 'earthMass' and 'planetEarthRads' and 'equilTemp' in row.keys():
        for p in referencePlanetSystems:
            similarityProf = p.similarityProfile
            xp = [p.mass, p.radius, p.eTemp, p.calculatedEqEarthTemp, p.calculatedEqVenusTemp, p.calculatedEqMarsTemp, p.calculatedEqTitanTemp, p.calculatedSurfaceEarthTemp, p.calculatedSurfaceVenusTemp,
                  p.calculatedSurfaceMarsTemp, p.calculatedSurfaceTitanTemp]
            x = [row['earthMass'], row['planetEarthRads'], row['equilTemp'], row['calculatedEqEarthTemp'], row['calculatedEqVenusTemp'], row['calculatedEqMarsTemp'], row['calculatedEqTitanTemp'],
                 row['calculatedSurfaceEarthTemp'],
                 row['calculatedSurfaceVenusTemp'], row['calculatedSurfaceMarsTemp'], row['calculatedSurfaceTitanTemp']]
            products = []
            for exponents in similarityProf.weightExponents:
                product = 1
                for i in range(len(exponents)):
                    product *= pow(1 - abs((xp[i] - x[i]) / (xp[i] + x[i])), exponents[i])
                    return product
            #     products.append(product)
            # return products
            #right now only looks at one set of exponents, could look at more

class planetSystem:
    def __init__(self, name, luminosity, distance, stellarRadius, albedo, beta, radius, mass, similarityProfile, emissivity=1):
        self.mass = mass  # in earth masses
        self.radius = radius  # earth radii
        self.name = name  # name of planet
        self.luminosity = luminosity  # in stellar units
        self.distance = distance  # in au
        self.stellarRadius = stellarRadius  # sun radii
        self.albedo = albedo
        self.beta = beta  # fraction of the planet that re-radiates the absorbed flux, default 1
        self.similarityProfile = similarityProfile
        self.emissivity = emissivity  # usually ~1
        self.eTemp = equilibriumTemp(self.luminosity, self.stellarRadius, self.distance, self.albedo, self.beta)
        self.calculatedEqEarthTemp = equilibriumTemp(self.luminosity, self.distance, self.stellarRadius, 0.301, 1)
        self.calculatedEqVenusTemp = equilibriumTemp(self.luminosity, self.distance, self.stellarRadius, 0.760, 1)
        self.calculatedEqMarsTemp = equilibriumTemp(self.luminosity, self.distance, self.stellarRadius, 0.250, 1)
        self.calculatedEqTitanTemp = equilibriumTemp(self.luminosity, self.distance, self.stellarRadius, 0.265, 1)
        self.calculatedSurfaceEarthTemp = surfaceTemp(0.385, self.luminosity, 1, self.distance, 0.301, 1)
        self.calculatedSurfaceVenusTemp = surfaceTemp(0.990, self.luminosity, 1, self.distance, 0.760, 1)
        self.calculatedSurfaceMarsTemp = surfaceTemp(0.073, self.luminosity, 1, self.distance, 0.250, 1)
        self.calculatedSurfaceTitanTemp = surfaceTemp(0.338, self.luminosity, 1, self.distance, 0.265, 1)


titan = planetSystem('Titan', 1, 9.5, 1, 0.265, 0.5, 0.404, 0.0225,
                     similarityProfile([[1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]])) #, [1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1], [1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]]))

# TITAN DECLARATION HERE D:




df = pd.read_csv('src\exoData.csv')


# provided data is as follows:
# Planet	T_(eq)	T_(s)	A_(b)	G_(n)	kappa
# ------------------------------------------------
# Venus	     229	730	    0.760	0.990	2.213
# Earth	     255	288	    0.301	0.385	1.033
# Mars	     210	214	    0.250	0.073	0.948
# Titan	     84	94	0.265	0.338	1.027

# outputted T_s as of 6/22/22:
# venus: 724.800
# earth: 287.558
# mars: 213.7265
# titan: 92.7514
# (pretty good)

# two different formulas for the same thing, for sanity check.




def calc():
    df['calculatedEqEarthTemp'] = df.apply(lambda row: equilibriumTemp(row['luminosity'], row['planetSemimajor'], row['stellarRadius'], 0.301, 1), axis=1)
    df['calculatedEqVenusTemp'] = df.apply(lambda row: equilibriumTemp(row['luminosity'], row['planetSemimajor'], row['stellarRadius'], 0.760, 1), axis=1)
    df['calculatedEqMarsTemp'] = df.apply(lambda row: equilibriumTemp(row['luminosity'], row['planetSemimajor'], row['stellarRadius'], 0.250, 1), axis=1)
    df['calculatedEqTitanTemp'] = df.apply(lambda row: equilibriumTemp(row['luminosity'], row['planetSemimajor'], row['stellarRadius'], 0.265, 1), axis=1)

    df['calculatedSurfaceEarthTemp'] = df.apply(lambda row: surfaceTemp(0.385, row['luminosity'], 1, row['planetSemimajor'], 0.301, 1), axis=1)
    df['calculatedSurfaceVenusTemp'] = df.apply(lambda row: surfaceTemp(0.990, row['luminosity'], 1, row['planetSemimajor'], 0.760, 1), axis=1)
    df['calculatedSurfaceMarsTemp'] = df.apply(lambda row: surfaceTemp(0.073, row['luminosity'], 1, row['planetSemimajor'], 0.250, 1), axis=1)
    df['calculatedSurfaceTitanTemp'] = df.apply(lambda row: surfaceTemp(0.338, row['luminosity'], 0.5, row['planetSemimajor'], 0.265, 1), axis=1)

    df['similarityToTitan'] = df.apply(lambda row: similarityIndex([titan], row),axis = 1)


print("old earth: " + str(surfaceTemp(0.385, 1, 1, 1, 0.301, 1)) + "\nold venus: " + str(surfaceTemp(0.990, 1, 1, 0.7233, 0.760, 1)) + "\nold mars: " + str(
    surfaceTemp(0.073, 1, 1, 1.5273, 0.250, 1)) + "\nold titan: " + str(surfaceTemp(0.338, 1, 1, 9.5, 0.265, 1)))
print("\nnew earth: " + str(newSurfaceTemp(0.385, 1, 1, 1, 0.301, 1)) + "\nnew venus: " + str(newSurfaceTemp(0.990, 1, 1, 0.7233, 0.760, 1)) + "\nnew mars: " + str(
    newSurfaceTemp(0.073, 1, 1, 1.5273, 0.250, 1)) + "\nnew titan: " + str(newSurfaceTemp(0.338, 1, 1, 9.5, 0.265, 1)))

calc()
# df.plot.scatter(x='planetEarthRads', y='similarityToTitan', c='equilTemp', colormap='viridis')
# plt.show()
# print(df.loc[0, :])

df.plot.scatter(x='planetEarthRads', y='similarityToTitan', c='equilTemp', colormap='viridis')
plt.show()
df.to_csv('src/out.csv')
display(df)
# we get some pretty big values here. sus
