import pandas as pd
from math import pi

pd.set_option('display.max_columns', None)

G = 6.6743 * pow(10, -11)
stephanBoltzmanConstant = 5.67 * pow(10, -8)

earthMass = 5.9722 * pow(10, 24)
earthRadius = 6371000
earthDensity = 5.51
earthEscapeVel = 11200
earthAlbedo = 0.3
earthGreenHouseConstant = 0.385
AU = 1.5 * pow(10, 8)

sunRadius = 696342


# Returns the density of the planet in earth densities given the mass in earth masses and the radius in earth radii
def calcDensity(mass, radius):
    volume = (4/3) * pi * pow(radius*earthRadius, 3)
    return ((mass*earthMass / volume) / 1000) / earthDensity


# Returns the escape velocity of the planet in earth escape velocities given the mass in earth masses and the radius in
# earth radii
def calcEscapeVel(eMass, eRadius):
    mass = eMass * earthMass
    radius = eRadius * earthRadius
    return pow((2 * G * mass) / radius, (1/2)) / earthEscapeVel


# Returns the surface temperature of the planet in degrees Kelvin given the radius of the star in km, the effective
# temperature of the star in degrees Kelvin, and the semi-major axis in km
def calcSurfTemp(starRadius, starEffTemp, semimajor):
    flux = (stephanBoltzmanConstant * pow(starRadius, 2) * pow(starEffTemp, 4)) / pow(semimajor, 2)
    print(flux)
    return pow(((1 - earthAlbedo) * flux) / ((4 * stephanBoltzmanConstant) * (1 - earthGreenHouseConstant)), (1/4))


# Returns the earth similarity index of a pair of planets given a list of the reference planet's attributes, a list of
# the comparative planet's corresponding attributes, and a list of the corresponding weights
def calcESI(refVals, vals, weights):
    esi = 1
    for i, val in enumerate(vals):
        esi *= pow(1 - abs((refVals[i] - val) / (refVals[i] + val)), weights[i])
    return esi


# Read in csv file
exoplanetData = pd.read_csv('exoData.csv')

# Drops the values I didn't need for this
exoplanetData = exoplanetData.drop('eccentricity', axis=1)
exoplanetData = exoplanetData.drop('earthFlux', axis=1)
exoplanetData = exoplanetData.drop('equilTemp', axis=1)
exoplanetData = exoplanetData.drop('stellarMass', axis=1)
exoplanetData = exoplanetData.drop('stLogSurfaceGrav', axis=1)
exoplanetData = exoplanetData.drop('Unnamed: 14', axis=1)
exoplanetData = exoplanetData.drop('luminosity', axis=1)

# Removes any exoplanets if its existence is controversial
hasData = exoplanetData['controversial'] == 0
exoplanetData = exoplanetData[hasData]

# Removes any exoplanets that don't have the required data
exoplanetData = exoplanetData[exoplanetData['planetSemimajor'].notna()]
exoplanetData = exoplanetData[exoplanetData['planetEarthRads'].notna()]
exoplanetData = exoplanetData[exoplanetData['earthMass'].notna()]
exoplanetData = exoplanetData[exoplanetData['stellarRadius'].notna()]
exoplanetData = exoplanetData[exoplanetData['stellarEffTemp'].notna()]

# Converts all stellar radii from solar radii to km and all planet semi-major axis from AU to km
exoplanetData['stellarRadius'] *= sunRadius
exoplanetData['planetSemimajor'] *= AU

# Adds and fills out the density, surface temperature, escape velocity, and earth similarity index columns
exoplanetData['density'] = calcDensity(exoplanetData['earthMass'], exoplanetData['planetEarthRads'])
exoplanetData['surfaceTemperature'] = calcSurfTemp(exoplanetData['stellarRadius'], exoplanetData['stellarEffTemp'],
                                                   exoplanetData['planetSemimajor'])
exoplanetData['escapeVelocity'] = calcEscapeVel(exoplanetData['earthMass'], exoplanetData['planetEarthRads'])
exoplanetData['ESI'] = calcESI([1, 1, 288, 1], [exoplanetData['planetEarthRads'], exoplanetData['escapeVelocity'],
                                                exoplanetData['surfaceTemperature'], exoplanetData['density']],
                               [1, 1, 1, 1])

# Prints out some relevant data
print(exoplanetData)
print(exoplanetData.loc[exoplanetData['surfaceTemperature'].idxmax()])
print(calcSurfTemp(sunRadius, 5772, AU))
