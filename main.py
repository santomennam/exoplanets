import pandas as pd
from math import pi
import os
from datetime import datetime
from matplotlib import pyplot as plt
from IPython.display import display
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

pd.set_option('display.max_columns', None)

# CONSTANTS
stefanBoltzmannConstant = 5.67e-8  # W m^-2 K^-4
earthEquilTemp = 278.5  # kelvin
AU = 1.5 * pow(10, 8)  # km
G = 6.67408 * pow(10, -11)  # m^3 kg^-1 s^-2
sunRadius = 696342  # km
solarLuminosity = 3.8282e26
earthMass = 5.9722 * pow(10, 24)  # kg
earthRadius = 6371000  # m
earthDensity = 5520  # kg m^-3
earthEscapeVel = 11200  # m s^-1


# Returns the density of the planet in earth densities, provided mass and radius are given in earth units
def calcDensity(mass: float, radius: float):
    volume = (4 / 3) * pi * pow(radius * earthRadius, 3)
    return (mass * earthMass / volume) / earthDensity


# calculate escape velocity in earth escape velocities, provided values in earth units
def calcEscapeVel(eMass: float, eRadius: float):
    # convert to kg
    mass = eMass * earthMass
    # convert to m
    radius = eRadius * earthRadius
    return pow((2 * G * mass) / radius, (1 / 2)) / earthEscapeVel


def calcSurfTemp(luminosity: float, semimajor: float, **kwargs):
    # first, we look for a solarBody in the kwargs. if one is provided, we use its albedo and greenhouse for the calculations
    if "solarBody" in kwargs.keys() is not None:
        p = kwargs["solarBody"]
        albedo, greenhouse, antiGreenhouse = p.albedo, p.greenhouse, p.antiGreenhouse
    # if no SolarBody was provided, we grab manually provided albedo and greenhouse values
    elif "albedo" in kwargs.keys() and "greenhouse" in kwargs.keys() and "antiGreenhouse" in kwargs.keys():
        albedo, greenhouse, antiGreenhouse = kwargs["albedo"], kwargs["greenhouse"], kwargs["antiGreenhouse"]
    else:
        raise Exception('Must provide a solar body or albedo, greenhouse, and antiGreenhouse')

    # solar luminosity will provided by exoplanet archive in units of log[solar luminosity] - we need to convert to watts
    luminosity = pow(10, luminosity) * solarLuminosity
    # semimajor is already in units of km, multiply by 1000 for m
    semimajor *= 10e3

    effectiveTemp = pow((1 - albedo) * luminosity / (16 * pi * semimajor ** 2 * stefanBoltzmannConstant), 1 / 4)

    surfaceTemp = effectiveTemp * pow(((1 - antiGreenhouse) * (1 + (3 / 4) * greenhouse) + antiGreenhouse / 2), 1 / 4)

    return surfaceTemp


# define a class for storing earth, venus, etc
class SolarBody:
    def __init__(self, name: str, albedo: float, greenhouseConstant: float, antiGreenhouse: float):
        self.name = name
        self.albedo = albedo
        self.greenhouse = greenhouseConstant
        self.antiGreenhouse = antiGreenhouse  # also known as grey opacity


# these represent our science fiction exoplanets
class CompBody:
    def __init__(self, name: str, semimajor, earthRads, escapeVel, surfaceTemp, density, exponents):
        self.semimajor = semimajor  # km
        self.radius = earthRads  # earth radii
        self.name = name  # name of planet
        self.escapeVel = escapeVel  # earth escape velocities
        self.surfTemp = surfaceTemp  # K
        self.density = density  # earth densities
        self.exponents = exponents  # pass in exponents for each CompBody that determine how planets will be compared to it
        self.refs = [self.radius, self.escapeVel, self.density, self.surfTemp, self.surfTemp, self.surfTemp,
                     self.surfTemp]


# perform the weighted product (ESI) calculation
def calcESI(refVals: list, vals: list, weights: list):
    esi = 1  # our ESI will start at 1 and will not get larger
    for i, val in enumerate(vals):
        esi *= pow(1 - abs((refVals[i] - val) / (refVals[i] + val)), weights[i])
    return round(esi, 3)


def pruneData(dataFrame: pd.DataFrame):
    """
    drop unnecessary columns and remove targets without enough information
    """

    # Removes any exoplanets if its existence is controversial
    dataFrame = dataFrame.loc[dataFrame['controversial'] == 0]
    # avoid multiple entries of the same planet by only using the best entry
    dataFrame = dataFrame.loc[dataFrame['defaultPlanetEntry'] == 1]

    # candidate not in soltype gets rid of kepler and tess candidates
    dataFrame = dataFrame[~dataFrame['solution'].str.contains('Candidate')]

    # remove eccentricities >= 1
    # dataFrame = dataFrame[~(dataFrame["eccentricity"] >= 0.1)]

    # Removes any exoplanets that don't have the required data
    dataFrame = dataFrame[dataFrame['planetSemimajor'].notna()]  # required for temperature calculation
    dataFrame = dataFrame[dataFrame['planetEarthRads'].notna()]  # required for density and escape velocity calculation
    dataFrame = dataFrame[dataFrame['earthMass'].notna()]  # required for density and escape velocity calculation
    dataFrame = dataFrame[dataFrame['stellarRadius'].notna()]  # required for temperature calculation
    dataFrame = dataFrame[dataFrame['stellarEffTemp'].notna()]  # required for temperature calculation
    dataFrame = dataFrame[dataFrame['luminosity'].notna()]  # required for temperature calculation

    return dataFrame


def dfCalcs(dataFrame):
    """
    Calculate density, surface temperatures, escape velocity
    :rtype: pd.DataFrame
    """
    # Converts all stellar radii from solar radii to km, and all planet semi-major axes from AU to km
    dataFrame['stellarRadius'] *= sunRadius
    dataFrame['planetSemimajor'] *= AU

    # Adds and fills out the density, surface temperature, and escape velocity columns
    dataFrame['density'] = calcDensity(dataFrame['earthMass'], dataFrame['planetEarthRads'])
    dataFrame['escapeVelocity'] = calcEscapeVel(dataFrame['earthMass'], dataFrame['planetEarthRads'])
    dataFrame['surfaceEarthTemp'] = dataFrame.apply(
        lambda row: calcSurfTemp(row['luminosity'], row['planetSemimajor'], solarBody=solarEarth), axis=1)
    dataFrame['surfaceVenusTemp'] = dataFrame.apply(
        lambda row: calcSurfTemp(row['luminosity'], row['planetSemimajor'], solarBody=solarVenus), axis=1)
    dataFrame['surfaceMarsTemp'] = dataFrame.apply(
        lambda row: calcSurfTemp(row['luminosity'], row['planetSemimajor'], solarBody=solarMars), axis=1)
    dataFrame['surfaceTitanTemp'] = dataFrame.apply(
        lambda row: calcSurfTemp(row['luminosity'], row['planetSemimajor'], solarBody=solarTitan), axis=1)

    return dataFrame


# calculate similarity to a given compbody for all exoplanets
def compareToPlanet(compBody: CompBody, dataframe, exponents: list):
    dataframe['similarityIndex_' + compBody.name] = dataframe.apply(
        lambda row: calcESI(compBody.refs, [row['planetEarthRads'],
                                            row['escapeVelocity'],
                                            row['density'],
                                            row['surfaceEarthTemp'],
                                            row['surfaceVenusTemp'],
                                            row['surfaceMarsTemp'],
                                            row['surfaceTitanTemp']],
                            exponents), axis=1)


def saveMatchesCSV(inFrame, value, path):
    """
    Saves the top three planets of a value 'value' and from a dataframe 'inFrame' to a csv in location 'path'
    :rtype: None
    """
    dataframe = inFrame.copy()
    dataframe.drop(dataframe.columns.difference(['Name', 'planetEarthRads', 'density', 'escapeVelocity',
                                                 'surfaceEarthTemp', value]), 1, inplace=True)
    outFrame = pd.DataFrame(columns=dataframe.columns)
    for i in range(3):
        highestIndex = dataframe[value].idxmax()
        tempDict = dataframe.loc[highestIndex].to_dict()
        outFrame = outFrame.append(tempDict, ignore_index=True)
        dataframe = dataframe.drop(highestIndex)
    outFrame.to_csv(path)


if __name__ == "__main__":
    # Creating SolarBodies. These will be our references for temperature calculations
    solarEarth = SolarBody('Earth', 0.301, 0.385, 0)
    solarVenus = SolarBody('Venus', 0.760, 0.990, 0)
    solarMars = SolarBody('Mars', 0.250, 0.073, 0)
    solarTitan = SolarBody('Titan', 0.265, 0.338, 4 / 7)

    # these constants will be used to for exponent normalization
    r = 1 / 3
    v = 1 / 3
    d = 1 / 3
    x = 1

    # Creating planets from real life or sci fi
    # name, semimajor, rads, escape, surface temp, density
    compEarth = CompBody('Earth', AU, 1, 1, 288, 1, [r, v, d, x, 0, 0, 0])
    compArrakis = CompBody('Arrakis', 2.3 * AU, 0.979, calcEscapeVel(0.843, 0.979), 325, calcDensity(0.843, 0.979),
                           [r, v, d, x, 0, 0, 0])
    compDhrawn = CompBody('Dhrawn', 0.3 * AU, 8.5, calcEscapeVel(3000, 8.5), 250, calcDensity(3000, 8.5),
                          [r, v, d, 0.5 * x, 0, 0.5 * x, 0])
    # TODO: placeholder exponents for Erid, need to discuss:
    compErid = CompBody('Erid', 0.224 * AU, 2.01, calcEscapeVel(8.47, 2.01), 480, calcDensity(8.47, 2.01),
                        [r, v, d, x, 0, 0, 0])
    compHekla = CompBody('Hekla', 0.3 * AU, 2, calcEscapeVel(1.45, 2), 260, calcDensity(1.45, 2), [r, v, d, 0, 0, x, 0])
    compMesklin = CompBody('Mesklin', 3.3 * AU, 7.85, calcEscapeVel(5086.51, 7.85), 100, calcDensity(5086.51, 7.85),
                           [r, v, d, 0, 0, 0, x])
    compSarr = CompBody('Sarr', 1.6 * AU, 0.76, calcEscapeVel(0.45, 0.76), 800, calcDensity(0.45, 0.76),
                        [r, v, d, 0, x, 0, 0])
    compTenebra = CompBody('Tenebra', 2 * AU, 3, calcEscapeVel(27, 3), 650, calcDensity(27, 3),
                           [r, v, d, 0.25 * x, 0.75 * x, 0, 0])

    # collect our solar bodies
    solarBodies = [solarEarth, solarVenus, solarMars, solarTitan]

    # collect our comparison (fictional+Earth) bodies
    compBodies = [compEarth, compArrakis, compErid, compMesklin, compDhrawn, compHekla, compSarr,
                  compTenebra]
    # collect our temperatures
    temps = [a.surfTemp for a in compBodies]

    # import and manage data
    df = pd.read_csv('2023.05.29_data_2.csv')
    df = pruneData(df)
    df = dfCalcs(df)

    # generate an output path, make directories
    path = 'outputs/' + str(datetime.strftime(datetime.now(), "%Y-%m-%d-%H-%M-%S")) + '/'
    os.mkdir(path)
    os.mkdir(path + '/plots')
    os.mkdir(path + '/csvs')

    # make a dataframe to store ESI values in
    compBodiesESI = pd.DataFrame(columns=["Name", "ESI"])

    # calculate the esi of each of the scifi bodies
    for compBody in compBodies:
        esi = calcESI(compEarth.refs, compBody.refs, compBody.exponents)
        row = {"Name": [compBody.name], "ESI": [esi]}
        compBodiesESI = compBodiesESI.append(row, ignore_index=True)
    display(compBodiesESI)
    compBodiesESI.to_csv(path + 'fictionalPlanetESI.csv')

    averageSimilarity = pd.DataFrame(columns=["Name", "Average"])
    # do the calculations
    for body in compBodies:
        name = body.name
        # make directories
        os.mkdir(path + '/plots/' + name)
        os.mkdir(path + '/csvs/' + name)
        # compare all the exoplanets to it
        compareToPlanet(body, df, body.exponents)
        # start collecting data with the "row" object
        row = {"Name": [body.name], "Average": [df['similarityIndex_' + name].mean()]}
        averageSimilarity = averageSimilarity.append(row, ignore_index=True)
        # plot and save figures
        plotPath = path + '/plots/' + name + '/'
        # radius vs similarity index
        df.plot.scatter(x='planetEarthRads', y=('similarityIndex_' + name), c='surfaceEarthTemp', colormap='plasma',
                        title=name)
        plt.savefig(plotPath + name + 'Radius.png')
        plt.close()
        # temperature vs similarity index
        df.plot.scatter(x='surfaceEarthTemp', y=('similarityIndex_' + name), c='planetEarthRads', colormap='viridis',
                        title=name)
        plt.savefig(plotPath + name + 'Temp.png')
        plt.close()
        # distance vs similarity index
        df.plot.scatter(x='distance', y=('similarityIndex_' + name), c='planetEarthRads', colormap='viridis',
                        title=name)
        plt.savefig(plotPath + name + 'Dist.png')
        plt.close()
        # orbital period vs similarity index
        df.plot.scatter(x='orbitalPeriod', y=('similarityIndex_' + name), c='planetEarthRads', colormap='viridis',
                        title=name)
        plt.savefig(plotPath + name + 'period.png')
        plt.close()
        # save to a csv
        df.to_csv(path + '/csvs/' + name + '/ref' + name + '.csv')
        # save top 5 matches to a csv
        saveMatchesCSV(df, 'similarityIndex_' + name, path + '/csvs/' + name + '/topFive' + name + '.csv')
    # save a master csv
    df.to_csv(path + "data" + '.csv')

    # sort the averageSimilarity dataframe by name and then save it
    averageSimilarity = averageSimilarity.sort_values('Name')
    averageSimilarity.to_csv(path + 'fictionalPlanetAverage.csv')
    print(averageSimilarity.to_string())
