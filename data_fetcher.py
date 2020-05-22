import requests
import numpy as np

national_data_source = "https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-json/dpc-covid19-ita-andamento-nazionale.json"
#Italy's COVID-19 repository with virus spread data, updated with daily resolution

def extract_sir(record, population: int = 0):
    '''
    Extracts SIR functions values given a record compliant to Italy COVID-19 repository json-based data format
    :param record: dict-like object (result of json conversion)
    :param population: total population = S(0) + index case
    :return: [susceptible,] infectious, recovered
    '''
    infectious = record.get("totale_positivi")
    recovered = record.get("dimessi_guariti") + record.get("deceduti")#- record.get("totale_positivi") - prev_r

    if population > 0:
        susceptible = population - infectious - recovered
        return susceptible, infectious, recovered
    else:
        return infectious, recovered


def fetch_data(src: str = national_data_source):
    '''
    Fetches data from Github repository. May throw a requests.exceptions.HTTPError
    :return: JSON decoded response only if HTTP request response-code is 200
    '''
    r = requests.get(src)
    r.raise_for_status()

    return r.json()


def process_data(data: np.ndarray, population: int = 0):
    '''
    Extracts array of SIR functions values off fetched data from Italy's COVID-19 opendata repository
    :param data: numpy array of given shape with fetched data from repository
    :param population: initial population or 1 + susceptible at time 0. default: 0, susceptible are not computed
    :return: numpy array of shape <data.shape> with SIR data
    '''
    if population > 0:
        sir = np.zeros(shape=(len(data), 3))
    else:
        sir = np.zeros(shape=(len(data), 2))

    for i, record in enumerate(data):
        sir[i] = extract_sir(record, population)

    return sir


def discretize(data: np.ndarray):
    '''
    Computes a discrete first derivative on extracted SIR functions values
    :param data: a numpy array of shape (x, 2|3)
    :return: array of differences in discrete domain
    '''
    diffs = np.diff(data, n=1, axis=0)

    return diffs


def plot(data: np.ndarray, discretized: bool = False):
    '''
    Plots data using matplotlib
    :param data: numpy array with SIR functions values
    :param discretized: if SIR functions are represented with derivates, integrate over values
    '''
    import matplotlib.pyplot as plt

    if discretized:
        plt.plot(data.cumsum(axis=0))
    else:
        plt.plot(data)

    if data.shape[1] == 3:
        plt.legend(('Susceptible', 'Infected', 'Recovered'))
    else:
        plt.legend(['Infected', 'Recovered'])

    plt.show()


def main():
    freshly_fetched = fetch_data()
    sir = process_data(freshly_fetched, population=250000)
    plot(sir)


if __name__ == '__main__':
    main()