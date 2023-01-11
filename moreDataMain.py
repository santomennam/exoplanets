import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import math
from math import pi
import os
from datetime import datetime
from matplotlib import pyplot as plt
from IPython.display import display

pd.set_option('display.max_columns', None)


#CONSTANTS
stefanBoltzmannConstant = 5.67e-8  # W m^-2 K^-4
earthEquilTemp = 278.5  # kelvin
AU = 1.5 * pow(10, 8)  # km
G = 6.67408 * pow(10, -11)  # m^3 kg^-1 s^-2

# DIMENSIONS
sunRadius = 696342  # km

earthMass = 5.9722 * pow(10, 24)  # kg
earthRadius = 6371000  # m
earthDensity = 5520  # kg m^-3
earthEscapeVel = 11200  # m s^-1


# Returns the density of the planet in earth densities given the mass in earth masses and the radius in earth radii
def calcDensity(mass, radius):
    volume = (4 / 3) * pi * pow(radius * earthRadius, 3)
    return (mass * earthMass / volume) / earthDensity


def calcEscapeVel(eMass, eRadius):
    mass = eMass * earthMass
    radius = eRadius * earthRadius
    return pow((2 * G * mass) / radius, (1 / 2)) / earthEscapeVel

    # see eqs 1 and 6 in Méndez, A. and Rivera-Valentín, E.G., 2017. The equilibrium temperature of planets in elliptical orbits. The Astrophysical Journal Letters, 837(1), p.L1.
def calcSurfTemp(stellarRadius, starEffTemp, semimajor, **kwargs):
    if "solarBody" in kwargs.keys() is not None:
        p = kwargs["solarBody"]
        container = (p.albedo, p.greenhouse)
    elif "albedo" in kwargs.keys() and "greenhouse" in kwargs.keys():
        container = (kwargs["albedo"], kwargs["greenhouse"])
    else:
        raise Exception('Must provide a solar body or albedo and greenhouse')
    flux = (stefanBoltzmannConstant * pow(stellarRadius, 2) * pow(starEffTemp, 4)) / pow(semimajor, 2)
    temp = pow(((1 - container[0]) * flux) / ((4 * stefanBoltzmannConstant) * (1 - container[1])), (1 / 4))
    return temp

def calcESI(refVals, vals, weights):
    esi = 1
    for i, val in enumerate(vals):
        esi *= pow(1 - abs((refVals[i] - val) / (refVals[i] + val)), weights[i])
    return esi


class solarBody:
    def __init__(self, name, albedo, greenhouseConstant):
        self.name = name
        self.albedo = albedo
        self.greenhouse = greenhouseConstant


class compBody:
    def __init__(self, name, semimajor, earthRads, escapeVel, surfaceTemp, density):
        self.semimajor = semimajor  # km
        self.radius = earthRads  # earth radii
        self.name = name  # name of planet
        self.escapeVel = escapeVel  # earth escape velocities
        self.surfTemp = surfaceTemp  # K
        self.density = density  # earth densities
        self.surfaceEarthTemp = 0  # K
        self.surfaceVenusTemp = 0  # K
        self.surfaceMarsTemp = 0  # K
        self.surfaceTitanTemp = 0  # K
        self.refs = []

    def surfaceTemps(self, earth, venus, mars, titan):
        self.surfaceEarthTemp = calcSurfTemp(sunRadius, 5870, self.semimajor, solarBody=earth)
        self.surfaceVenusTemp = calcSurfTemp(sunRadius, 5870, self.semimajor, solarBody=venus)
        self.surfaceMarsTemp = calcSurfTemp(sunRadius, 5870, self.semimajor, solarBody=mars)
        self.surfaceTitanTemp = calcSurfTemp(sunRadius, 5870, self.semimajor, solarBody=titan)
        self.refs = [self.radius, self.escapeVel, self.density, self.surfaceEarthTemp, self.surfaceVenusTemp,
                     self.surfaceMarsTemp, self.surfaceTitanTemp]

    def setSurfaceTemp(self, temp):
        self.surfaceEarthTemp = self.surfaceVenusTemp = self.surfaceMarsTemp = self.surfaceTitanTemp = temp
        self.refs = [self.radius, self.escapeVel, self.density, self.surfaceEarthTemp, self.surfaceVenusTemp,
                     self.surfaceMarsTemp, self.surfaceTitanTemp]

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
    print("Rows before prune: ", str(df.shape[0]))
    # Removes any exoplanets that don't have the required data
    df = df[df['planetSemimajor'].notna()]
    df = df[df['planetEarthRads'].notna()]
    df = df[df['earthMass'].notna()]
    df = df[df['stellarRadius'].notna()]
    df = df[df['stellarEffTemp'].notna()]

    df = df[df['controversial'] == 0]
    df = df[~df['solution'].str.contains('Candidate')] #candidate not in soltype gets rid of kepler and tess candidates

    print("Rows after prune: ", str(df.shape[0]))
    return df


def dfCalcs(dataFrame):
    # Converts all stellar radii from solar radii to km and all planet semi-major axis from AU to km
    dataFrame['stellarRadius'] *= sunRadius
    dataFrame['planetSemimajor'] *= AU

    # Adds and fills out the density, surface temperature, and escape velocity columns
    dataFrame['density'] = calcDensity(dataFrame['earthMass'], dataFrame['planetEarthRads'])
    dataFrame['escapeVelocity'] = calcEscapeVel(dataFrame['earthMass'], dataFrame['planetEarthRads'])
    dataFrame['surfaceEarthTemp'] = dataFrame.apply(lambda row: calcSurfTemp(row['stellarRadius'], row['stellarEffTemp'], row['planetSemimajor'], solarBody=solarEarth), axis=1)
    dataFrame['surfaceVenusTemp'] = dataFrame.apply(lambda row: calcSurfTemp(row['stellarRadius'], row['stellarEffTemp'], row['planetSemimajor'], solarBody=solarVenus), axis=1)
    dataFrame['surfaceMarsTemp'] = dataFrame.apply(lambda row: calcSurfTemp(row['stellarRadius'], row['stellarEffTemp'], row['planetSemimajor'], solarBody=solarMars), axis=1)
    dataFrame['surfaceTitanTemp'] = dataFrame.apply(lambda row: calcSurfTemp(row['stellarRadius'], row['stellarEffTemp'], row['planetSemimajor'], solarBody=solarTitan), axis=1)

    return dataFrame


def compareToPlanet(compBody, dataframe, exponents):
    df = dataframe
    # print("solarBody refs: ", solarBody.refs)
    df['similarityIndex_' + compBody.name] = df.apply(lambda row: calcESI(compBody.refs, [row['planetEarthRads'], row['escapeVelocity'], row['density'], row['surfaceEarthTemp'], row['surfaceVenusTemp'], row['surfaceMarsTemp'], row['surfaceTitanTemp']], exponents), axis=1)

# Dumps the top five planets of a value 'value' and from a dataframe 'inFrame' to a csv in location 'path'
def getMatchesCSV(inFrame, value, path):
    dataframe = inFrame.copy()
    dataframe.drop(dataframe.columns.difference(['Name', 'planetEarthRads', 'density', 'escapeVelocity',
                                                 'surfaceEarthTemp', value]), 1, inplace=True)
    outFrame = pd.DataFrame(columns=dataframe.columns)
    for i in range(5):
        highestIndex = dataframe[value].idxmax()
        tempDict = dataframe.loc[highestIndex].to_dict()
        outFrame = outFrame.append(tempDict, ignore_index=True)
        dataframe = dataframe.drop(highestIndex)
    outFrame.to_csv(path)


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
compErid = compBody('Erid', 0.224*AU, 2.01, calcEscapeVel(8.47, 2.01), 480, calcDensity(8.47, 2.01))
compMesklin = compBody('Mesklin', 3.3*AU, 7.85, calcEscapeVel(5086.51, 7.85), 100, calcDensity(5086.51, 7.85))
compHabranah = compBody('Habranah', 1.5*AU, 0.21, calcEscapeVel(0.01, 0.21), 195, calcDensity(0.01, 0.21))
compDhrawn = compBody('Dhrawn', 0.3*AU, 8.5, calcEscapeVel(3000, 8.5), 250, calcDensity(3000, 8.5))
compHekla = compBody('Hekla', 0.3*AU, 2, calcEscapeVel(1.45, 2), 260, calcDensity(1.45, 2))
compSarr = compBody('Sarr', 1.6*AU, 0.76, calcEscapeVel(0.45, 0.76), 800, calcDensity(0.45, 0.76))
compTenebra = compBody('Tenebra', 2*AU, 3, calcEscapeVel(27, 3), 650, calcDensity(27, 3))

solarBodies = [solarEarth, solarVenus, solarMars, solarTitan]

compBodies = [compEarth, compJupiter, compArrakis, compErid, compMesklin, compHabranah, compDhrawn, compHekla, compSarr,
              compTenebra]  # , Arrakis]



for body in compBodies:
    body.setSurfaceTemp(body.surfTemp)  # MUST CALL! -> this can be expedited, is an explicit call in case of surface temperature calculations being performed on the fly

# import and manage data
df = pd.read_csv('evenMoreExoData.csv')
df = pruneData(df)
df = dfCalcs(df)

exponents = [1, 1, 1, 1, 0, 0, 0]
path = 'outputs/' + "more "+str(datetime.strftime(datetime.now(), "%Y-%m-%d-%H-%M")) + '/'
os.mkdir(path)
os.mkdir(path + '/plots')
os.mkdir(path + '/csvs')
compBodiesESI = pd.DataFrame(columns=["Name","ESI"])
for compBody in compBodies:
    esi = calcESI(compEarth.refs, compBody.refs, exponents)
    row = {"Name": [compBody.name], "ESI": [esi]}
    compBodiesESI = compBodiesESI.append(row, ignore_index=True)
display(compBodiesESI)
compBodiesESI.to_csv(path + 'fictionalPlanetESI.csv')

averageSimilarity = pd.DataFrame(columns=["Name","Average"])
for body in compBodies:
    name = body.name
    os.mkdir(path + '/plots/' + name)
    os.mkdir(path + '/csvs/' + name)

    compareToPlanet(body, df, exponents)
    plotPath = path + '/plots/' + name + '/'
    df.plot.scatter(x='planetEarthRads', y=('similarityIndex_' + name),  c='surfaceEarthTemp', colormap='plasma', title=name)
    plt.savefig(plotPath + name + 'Radius.png')
    plt.close()
    df.plot.scatter(x='surfaceEarthTemp', y=('similarityIndex_' + name), c='planetEarthRads', colormap='viridis', title=name)
    plt.xlim(0, 2800) #  these are just manually discarding outliers for aesthetic reasons
    plt.savefig(plotPath + name + 'Temp.png')
    plt.close()
    df.plot.scatter(x='distance', y=('similarityIndex_' + name), c='planetEarthRads', colormap='viridis', title=name)
    plt.xlim(0, 3000)
    plt.savefig(plotPath + name + 'Dist.png')
    plt.close()
    df.plot.scatter(x='orbitalPeriod', y=('similarityIndex_' + name), c='planetEarthRads', colormap='viridis', title=name)
    plt.xlim(0, 400)
    plt.savefig(plotPath + name + 'period.png')
    plt.close()
    row = {"Name": [body.name], "Average": [df['similarityIndex_' + name].mean()]}
    averageSimilarity = averageSimilarity.append(row, ignore_index=True)
    getMatchesCSV(df, 'similarityIndex_' + name, path + '/csvs/' + name + '/topFive' + name + '.csv')
df.to_csv(path + "data" + '.csv')
averageSimilarity = averageSimilarity.sort_values('Name')
averageSimilarity.to_csv(path + 'fictionalPlanetAverage.csv')

