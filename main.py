import pandas as pd
import math
from math import pi
import os
from datetime import datetime
from matplotlib import pyplot as plt
from IPython.display import display

pd.set_option('display.max_columns', None)

stefanBoltzmannConstant = 5.67e-8
earthDistFromSun = 149.6e9  # km
earthEquilTemp = 278.5  # kelvin
AU = 1.5 * pow(10, 8)
G = 6.67408 * pow(10, -11)

# CONSTANTS
sunRadius = 696342

earthMass = 5.9722 * pow(10, 24)
earthRadius = 6371000
earthDensity = 5.515
earthEscapeVel = 11200



# Returns the density of the planet in earth densities given the mass in earth masses and the radius in earth radii
def calcDensity(mass, radius):
    volume = (4 / 3) * pi * pow(radius * earthRadius, 3)
    return ((mass * earthMass / volume) / 1000) / earthDensity


def calcEscapeVel(eMass, eRadius):
    mass = eMass * earthMass
    radius = eRadius * earthRadius
    return pow((2 * G * mass) / radius, (1 / 2)) / earthEscapeVel

    # see eqs 1 and 2 in Méndez, A. and Rivera-Valentín, E.G., 2017. The equilibrium temperature of planets in elliptical orbits. The Astrophysical Journal Letters, 837(1), p.L1.


def calcSurfTemp(stellarRadius, starEffTemp, semimajor, **kwargs):
    # print("in calculate surface temp. stellarRadius: ", stellarRadius, " starEffTemp: ", starEffTemp, " semimajor: ", semimajor)
    if "solarBody" in kwargs.keys() is not None:
        p = kwargs["solarBody"]
        container = (p.albedo, p.greenhouse)
    elif "albedo" in kwargs.keys() and "greenhouse" in kwargs.keys():
        container = (kwargs["albedo"], kwargs["greenhouse"])
    else:
        raise Exception('Must provide a solar body or albedo and greenhouse')
    # print("container: ", container)
    flux = (stefanBoltzmannConstant * pow(stellarRadius, 2) * pow(starEffTemp, 4)) / pow(semimajor, 2)
    temp = pow(((1 - container[0]) * flux) / ((4 * stefanBoltzmannConstant) * (1 - container[1])), (1 / 4))
    # print("temp: ", temp)
    return temp

# list of exponents for: mass, radius, eTemp, calculatedEqEarthTemp, calculatedEqVenusTemp, calculatedEqMarsTemp, calculatedEqTitanTemp, calculatedSurfaceEarthTemp,
# calculatedSurfaceVenusTemp, calculatedSurfaceMarsTemp, calculatedSurfaceTitanTemp

# sample exponents:
# earth:         [1,1,1,1,0,0,0,1,0,0,0]
# venus:         [1,1,1,0,1,0,0,0,1,0,0]
# mars:          [1,1,1,0,0,1,0,0,0,1,0]
# titan:         [1,1,1,0,0,0,1,0,0,0,1]
# no temp:       [1,1,1,0,0,0,0,0,0,0,0]
# body-neutral:  [1,1,1,1,1,1,1,1,1,1,1]


#  SI = product(i=1->m) of (1-abs((xp_i-x_i)/(xp_i+x_i)))^w_i
# where xp is the planetary reference value, x_i is the actual value of the comparison planet, and m is the number of properties.
# takes a list of planets, each with a list of list of exponents and their own characteristics. these are the base planets.
# also takes a dataframe row and compares the reference planets to the real database planet.


def calcESI(refVals, vals, weights):
    esi = 1
    # print("refVals: ", refVals, " vals: ", vals, " weights: ", weights, "\n")
    for i, val in enumerate(vals):
        esi *= pow(1 - abs((refVals[i] - val) / (refVals[i] + val)), weights[i])
    # # print(refVals, " , ", vals, ", ", weights, ", ", esi)
    # print("esi: ", esi)
    return esi


class solarBody:
    def __init__(self, name, albedo, greenhouseConstant):
        self.name = name
        self.albedo = albedo
        self.greenhouse = greenhouseConstant



class compBody:
    def __init__(self, name, semimajor, earthRads, escapeVel, surfaceTemp, density):
        self.semimajor = semimajor #in km
        self.radius = earthRads  # earth radii
        self.name = name  # name of planet
        self.escapeVel = escapeVel  # escape velocity in earth escape vels
        self.surfTemp = surfaceTemp  # surface temperature in kelvin
        self.density = density  # density in earth densities

    def surfaceTemps(self, earth, venus, mars, titan):
        self.surfaceEarthTemp = calcSurfTemp(sunRadius, 5870, self.semimajor, solarBody=earth)
        self.surfaceVenusTemp = calcSurfTemp(sunRadius, 5870, self.semimajor, solarBody=venus)
        self.surfaceMarsTemp = calcSurfTemp(sunRadius, 5870, self.semimajor, solarBody=mars)
        self.surfaceTitanTemp = calcSurfTemp(sunRadius, 5870, self.semimajor, solarBody=titan)
        self.refs = [self.radius, self.escapeVel, self.density, self.surfaceEarthTemp, self.surfaceVenusTemp, self.surfaceMarsTemp, self.surfaceTitanTemp]

    def setSurfaceTemp(self, temp):  # this is gimmicky as hell
        self.eTemp, self.calculatedEqEarthTemp, self.calculatedEqVenusTemp, self.calculatedEqMarsTemp, self.calculatedEqTitanTemp, self.calculatedSurfaceEarthTemp, self.calculatedSurfaceVenusTemp, self.calculatedSurfaceMarsTemp, self.calculatedSurfaceTitanTemp = [
                                                                                                                                                                                                                                                                           temp] * 9

# provided data is as follows:
# Planet    T_(eq)    T_(s)    A_(b)    G_(n)    kappa
# ------------------------------------------------
# Venus         229    730        0.760    0.990    2.213
# Earth         255    288        0.301    0.385    1.033
# Mars         210    214        0.250    0.073    0.948
# Titan         84    94    0.265    0.338    1.027

# outputted T_s as of 6/22/22:
# venus: 724.800
# earth: 287.558
# mars: 213.7265
# titan: 92.7514
# (pretty good)

def pruneData(df):
    df = df.drop('eccentricity', axis=1)
    df = df.drop('earthFlux', axis=1)
    df = df.drop('equilTemp', axis=1)
    df = df.drop('stellarMass', axis=1)
    df = df.drop('stLogSurfaceGrav', axis=1)
    df = df.drop('Unnamed: 14', axis=1)
    df = df.drop('luminosity', axis=1)

    # Removes any exoplanets if its existence is controversial
    hasData = df['controversial'] == 0
    df = df[hasData]

    # Removes any exoplanets that don't have the required data
    df = df[df['planetSemimajor'].notna()]
    df = df[df['planetEarthRads'].notna()]
    df = df[df['earthMass'].notna()]
    df = df[df['stellarRadius'].notna()]
    df = df[df['stellarEffTemp'].notna()]

    return df


def dfCalcs(df):
    # Converts all stellar radii from solar radii to km and all planet semi-major axis from AU to km
    df['stellarRadius'] *= sunRadius
    df['planetSemimajor'] *= AU

    # Adds and fills out the density, surface temperature, escape velocity, and earth similarity index columns
    df['density'] = calcDensity(df['earthMass'], df['planetEarthRads'])
    df['escapeVelocity'] = calcEscapeVel(df['earthMass'], df['planetEarthRads'])
    df['surfaceEarthTemp'] = df.apply(lambda row: calcSurfTemp(row['stellarRadius'], row['stellarEffTemp'], row['planetSemimajor'], solarBody=solarEarth), axis=1)
    df['surfaceVenusTemp'] = df.apply(lambda row: calcSurfTemp(row['stellarRadius'], row['stellarEffTemp'], row['planetSemimajor'], solarBody=solarVenus), axis=1)
    df['surfaceMarsTemp'] = df.apply(lambda row: calcSurfTemp(row['stellarRadius'], row['stellarEffTemp'], row['planetSemimajor'], solarBody=solarMars), axis=1)
    df['surfaceTitanTemp'] = df.apply(lambda row: calcSurfTemp(row['stellarRadius'], row['stellarEffTemp'], row['planetSemimajor'], solarBody=solarTitan), axis=1)

    return df


############ NEEDS WORK FROM HERE DOWN ##############


def compareToPlanet(compBody, dataframe, exponents):
        df = dataframe
        # print("solarBody refs: ", solarBody.refs)
        df['similarityIndex_' + compBody.name] = df.apply(lambda row: calcESI(compBody.refs, [row['planetEarthRads'], row['escapeVelocity'], row['density'], row['surfaceEarthTemp'], row['surfaceVenusTemp'], row['surfaceMarsTemp'], row['surfaceTitanTemp']], exponents), axis=1)

# when creating planets, enter 0 for stellar radius, distance, luminosity, and albedo if unknown.
# if known, call the "calcs" function on the planet object after creation
# Creating solarBodies for references
solarEarth = solarBody('Earth', 0.301, 0.385)
solarVenus = solarBody('Venus', 0.760, 0.990)
solarMars = solarBody('Mars', 0.250, 0.073)
solarTitan = solarBody('Titan', 0.265, 0.338)

# Creating planets from real life or sci fi
# name, semimajor, rads, escape, surface temp, density
compEarth = compBody('Earth', AU, 1, 1, 288, 1)
compJupiter = compBody('Jupiter', 5.2038*AU, 11.2, 5.317247, 123, 0.2404)
compArrakis = compBody('Arrakis', 2.3*AU, 0.979, calcEscapeVel(0.843, 0.979), 325, calcDensity(0.843, 0.979))
compErid = compBody('Erid', 0.4*AU, 2, calcEscapeVel(8, 2), 420, calcDensity(8, 2))
compMesklin = compBody('Mesklin', 3.3*AU, 7.85, calcEscapeVel(5086.51, 7.85), 100, calcDensity(5086.51, 7.85))
compHabranah = compBody('Habranah', 1.5*AU, 0.21, calcEscapeVel(0.01, 0.21), 195, calcDensity(0.01, 0.21))
compDhrawn = compBody('Dhrawn', 0.3*AU, 8.5, calcEscapeVel(3000, 8.5), 250, calcDensity(3000, 8.5))
compHekla = compBody('Hekla', 0.3*AU, 2, calcEscapeVel(1.45, 2), 260, calcDensity(1.45, 2))
compSarr = compBody('Sarr', 1.6*AU, 0.76, calcEscapeVel(0.45, 0.76), 800, calcDensity(0.45, 0.76))
compTenebra = compBody('Tenebra', 2*AU, 3, calcEscapeVel(27, 3), 650, calcDensity(27, 3))


# Arrakis = planetSystem('Arrakis', 0, 0, 0, 0, 1, 0.2723366, 0.0667)
# Arrakis.setSurfaceTemp(325)  # from the wiki or literature # star is canopus

# collecting

solarBodies = [solarEarth, solarVenus, solarMars, solarTitan]


compBodies = [compEarth, compJupiter, compArrakis, compErid, compMesklin, compHabranah, compDhrawn, compHekla, compSarr, compTenebra]  # , Arrakis]


for body in compBodies:
    body.surfaceTemps(solarEarth, solarVenus, solarMars, solarTitan) #MUST CALL!
print(compBodies[0].surfaceEarthTemp)
# import and manage data
df = pd.read_csv('exoData.csv')
df = pruneData(df)
df = dfCalcs(df)
print(df)


exponents = [1, 1, 1, 1, 0, 0, 0]
path = 'outputs/' + str(datetime.strftime(datetime.now(), "%Y-%m-%d-%H-%M")) + '/'
os.mkdir(path)
os.mkdir(path + '/plots')
os.mkdir(path + '/csvs')

for body in compBodies:
    name = body.name
    os.mkdir(path + '/plots/' + name)
    os.mkdir(path + '/csvs/' + name)
    compareToPlanet(body, df, exponents)
    df.plot.scatter(x='planetEarthRads', y=('similarityIndex_' + name), title=name)
    plt.savefig(path + '/plots/' + name + '/ref' + name + '.png')
    df.to_csv(path + '/csvs/' + name + '/ref' + name + '.csv')

# display(df)
# we get some pretty big values here. sus
