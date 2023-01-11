import pandas as pd
from math import pi
import os
from datetime import datetime
from matplotlib import pyplot as plt
from IPython.display import display
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

pd.set_option('display.max_columns', None)

stefanBoltzmannConstant = 5.67e-8  # W m^-2 K^-4
earthEquilTemp = 278.5  # kelvin
AU = 1.5 * pow(10, 8)  # km
G = 6.67408 * pow(10, -11)  # m^3 kg^-1 s^-2

# CONSTANTS
sunRadius = 696342  # km

earthMass = 5.9722 * pow(10, 24)  # kg
earthRadius = 6371000  # m
earthDensity = 5520  # kg m^-3
earthEscapeVel = 11200  # m s^-1

# Returns the density of the planet in earth densities, provided mass and radius are given in earth units
def calcDensity(mass, radius):
    volume = (4 / 3) * pi * pow(radius * earthRadius, 3)
    return (mass * earthMass / volume) / earthDensity

# calc escape vel provided values in earth units
def calcEscapeVel(eMass, eRadius):
    mass = eMass * earthMass
    radius = eRadius * earthRadius
    return pow((2 * G * mass) / radius, (1 / 2)) / earthEscapeVel


def calcSurfTemp(stellarRadius, starEffTemp, semimajor, **kwargs):
    # first , we look for a solarBody in the kwargs. if one is provided, we use its albedo and greenhouse for the calculations
    if "solarBody" in kwargs.keys() is not None:
        p = kwargs["solarBody"]
        albedo, greenhouse = p.albedo, p.greenhouse
    #here, we grab manually provided albedo and greenhouse if no solarBody was provided
    elif "albedo" in kwargs.keys() and "greenhouse" in kwargs.keys():
        albedo, greenhouse = kwargs["albedo"], kwargs["greenhouse"]
    #uh oh
    else:
        raise Exception('Must provide a solar body or albedo and greenhouse')
    #this calculates incident stellar flux from eq 1 in mendez et al
    flux = (stefanBoltzmannConstant * pow(stellarRadius, 2) * pow(starEffTemp, 4)) / pow(semimajor, 2)
    #slightly modified version of eqs 2, 4, 5, and 6
    temp = pow(((1 - albedo) * flux) / ((4 * stefanBoltzmannConstant) * (1 - greenhouse)), (1 / 4))
    # print("temp: ", temp)
    return temp

#perform the weighted product (ESI) calc
def calcESI(refVals, vals, weights):
    esi = 1
    for i, val in enumerate(vals):
        esi *= pow(1 - abs((refVals[i] - val) / (refVals[i] + val)), weights[i])
    return esi

#define a class for storing earth, venus, etc
class solarBody:
    def __init__(self, name, albedo, greenhouseConstant):
        self.name = name
        self.albedo = albedo
        self.greenhouse = greenhouseConstant

#these represent our science fiction exoplanets
class compBody:
    def __init__(self, name, semimajor, earthRads, escapeVel, surfaceTemp, density, exponents):
        self.semimajor = semimajor  # km
        self.radius = earthRads  # earth radii
        self.name = name  # name of planet
        self.escapeVel = escapeVel  # earth escape velocities
        self.surfTemp = surfaceTemp  # K
        self.density = density  # earth densities
        self.exponents = exponents #pass in exponents for each compBody that determine how planets will be compared to it
        self.refs = [self.radius, self.escapeVel, self.density, self.surfTemp, self.surfTemp, self.surfTemp, self.surfTemp] #ignore how horribly inelegant this is


# provided data is as follows:
# Planet    T_(eq)    T_(s)    A_(b)    G_(n)    kappa
# ------------------------------------------------
# Venus         229    730        0.760    0.990    2.213
# Earth         255    288        0.301    0.385    1.033
# Mars         210    214        0.250    0.073    0.948
# Titan         84    94    0.265    0.338    1.027


def pruneData(dataFrame):
    #drop unnecessary columns
    dataFrame = dataFrame.drop('eccentricity', axis=1)
    dataFrame = dataFrame.drop('earthFlux', axis=1)
    dataFrame = dataFrame.drop('equilTemp', axis=1)
    dataFrame = dataFrame.drop('stellarMass', axis=1)
    dataFrame = dataFrame.drop('stLogSurfaceGrav', axis=1)
    dataFrame = dataFrame.drop('Unnamed: 14', axis=1)
    dataFrame = dataFrame.drop('luminosity', axis=1)

    # Removes any exoplanets if its existence is controversial
    hasData = dataFrame['controversial'] == 0
    dataFrame = dataFrame[hasData]

    # Removes any exoplanets that don't have the required data
    dataFrame = dataFrame[dataFrame['planetSemimajor'].notna()]
    dataFrame = dataFrame[dataFrame['planetEarthRads'].notna()]
    dataFrame = dataFrame[dataFrame['earthMass'].notna()]
    dataFrame = dataFrame[dataFrame['stellarRadius'].notna()]
    dataFrame = dataFrame[dataFrame['stellarEffTemp'].notna()]

    return dataFrame


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


#calculate similarity to a given compbody for all exoplanets
def compareToPlanet(compBody, dataframe, exponents):
    df = dataframe
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
solarErid = solarBody('Erid',0.014,0)

# Creating planets from real life or sci fi
# name, semimajor, rads, escape, surface temp, density
compEarth = compBody('Earth', AU, 1, 1, 288, 1)
# compJupiter = compBody('Jupiter', 5.2038*AU, 11.2, 5.317247, 123, 0.2404)
compArrakis = compBody('Arrakis', 2.3*AU, 0.979, calcEscapeVel(0.843, 0.979), 325, calcDensity(0.843, 0.979))
compDhrawn = compBody('Dhrawn', 0.3*AU, 8.5, calcEscapeVel(3000, 8.5), 250, calcDensity(3000, 8.5))
compErid = compBody('Erid', 0.224*AU, 2.01, calcEscapeVel(8.47, 2.01), 480, calcDensity(8.47, 2.01))
# compHabranah = compBody('Habranah', 1.5*AU, 0.21, calcEscapeVel(0.01, 0.21), 195, calcDensity(0.01, 0.21))
compHekla = compBody('Hekla', 0.3*AU, 2, calcEscapeVel(1.45, 2), 260, calcDensity(1.45, 2))
compMesklin = compBody('Mesklin', 3.3*AU, 7.85, calcEscapeVel(5086.51, 7.85), 100, calcDensity(5086.51, 7.85))
compSarr = compBody('Sarr', 1.6*AU, 0.76, calcEscapeVel(0.45, 0.76), 800, calcDensity(0.45, 0.76))
compTenebra = compBody('Tenebra', 2*AU, 3, calcEscapeVel(27, 3), 650, calcDensity(27, 3))

#collect our solar bodies
solarBodies = [solarEarth, solarVenus, solarMars, solarTitan]
#collect our comparison (scifi) bodies
compBodies = [compEarth, compArrakis, compErid, compMesklin, compDhrawn, compHekla, compSarr,
              compTenebra] #compHabranah and compJupiter eliminated 1/11/23 by ben and sage

# import and manage data
df = pd.read_csv('exoData.csv')
df = pruneData(df)
df = dfCalcs(df)

#obviously these need tuning

#generate an output path, make directories
path = 'outputs/' + str(datetime.strftime(datetime.now(), "%Y-%m-%d-%H-%M-%S")) + '/'
os.mkdir(path)
os.mkdir(path + '/plots')
os.mkdir(path + '/csvs')

#make a dataframe to store ESI values in
compBodiesESI = pd.DataFrame(columns=["Name","ESI"])

#calculate the esi of each of the scifi bodies
for compBody in compBodies:
    esi = calcESI(compEarth.refs, compBody.refs, compBody.exponents)
    row = {"Name": [compBody.name], "ESI": [esi]}
    compBodiesESI = compBodiesESI.append(row, ignore_index=True)
display(compBodiesESI)
compBodiesESI.to_csv(path + 'fictionalPlanetESI.csv')

averageSimilarity = pd.DataFrame(columns=["Name","Average"])
#do the calculations
for body in compBodies:
    # **** for each of our sci fi planets, we: ----------------------------------
    name = body.name
    # **** make directories  ----------------------------------------------------
    os.mkdir(path + '/plots/' + name)
    os.mkdir(path + '/csvs/' + name)
    # **** compare all the exoplanets to it -------------------------------------
    compareToPlanet(body, df, body.exponents)
        # stats = df['similarityIndex_' + name].describe()
        # print(stats)
    # **** start collecting data with the "row" object --------------------------
    row = {"Name": [body.name], "Average": [df['similarityIndex_' + name].mean()]}
    averageSimilarity = averageSimilarity.append(row, ignore_index=True)
    # **** plot and save figs ---------------------------------------------------
    df.plot.scatter(x='surfaceEarthTemp', y=('similarityIndex_' + name),  c='planetEarthRads', colormap='viridis', title=name)
    plt.savefig(path + '/plots/' + name + '/ref' + name + '.png')
    # **** save to a csv --------------------------------------------------------
    df.to_csv(path + '/csvs/' + name + '/ref' + name + '.csv')
    # **** save top 5 matches to a csv ------------------------------------------
    getMatchesCSV(df, 'similarityIndex_' + name, path + '/csvs/' + name + '/topFive' + name + '.csv')

# **** sort the averageSimilarity by name and then save it ----------------------
averageSimilarity = averageSimilarity.sort_values('Name')
averageSimilarity.to_csv(path + 'fictionalPlanetAverage.csv')
