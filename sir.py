import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize


class SIRPredict:
    '''
    Represents the SIR model facilities to compute and predict the spread of
    a virus given the value of i
    '''

    def __init__(self, population: int, beta: float, gamma: float, t: int = 0, initsir: tuple = None):
        self.beta = beta
        self.gamma = gamma

        self.t = 0                              #Initial time
        self.n = population

        if not initsir is None:
            self.s, self.i, self.r = initsir
        else:
            self.reset_model()

    def reset_model(self):
        '''
        Initialize and resets model to t=0.
        Susceptible = population - 1
        Infected = 1 (index case)
        '''
        self.s = self.n - 1
        self.i = 1
        self.r = 0

    def solve(self):
        dsdt = -(self.beta * self.s * self.i / self.n)
        didt = (self.beta * self.s * self.i / self.n) - self.gamma * self.i
        drdt = self.gamma * self.i

        self.s += dsdt
        self.i += didt
        self.r += drdt

        return self.s, self.i, self.r

    def predict(self, day: int):
        '''
        Predict epidemic spread given the day, in the discrete domain.
        :param day: Day from the beginning of spread
        :return: Tuple with SIR values in day t
        '''
        for _ in range(day - self.t):
            self.solve()

        return self.s, self.i, self.r

    def spread_predict(self, final_day: int):
        '''
        Predict epidemic spread for a timespan, in the discrete domain.
        :param day: Final day of prediction
        :return: Numpy array with SIR values of predicted spread
        '''

        predicted = np.ndarray(shape=(final_day, 3))
        for i in range(final_day):
            predicted[i] = self.solve()

        return predicted

class SIRInterpolation:
    def __init__(self, sir_values: np.ndarray):
        '''
        Initializes the SIR model given an array of SIR functions first discrete derivative.
        :param sir_values: Array with SIR values.
        '''
        self.data = sir_values
        self.n = sir_values[0][0] + 1

    def loss_rms(self, point, data):
        '''
        Tries given
        :param point:
        :param data:
        :return: float
        '''

        size = len(data)
        beta, gamma = point

        def next_dt_SIR(t, sir_data):
            s, i, r = sir_data
            return (-beta*s*i/self.n, beta*s*i/self.n-gamma*i, gamma*i)

        solution = solve_ivp(next_dt_SIR,
                             (0, size),                     #Integration interval
                             self.data[0],                  #Initial state, float array of 3
                             t_eval=np.arange(0, size, 1),  #Discrete time interval
                             vectorized=True                #Functions supports vectors
                             )
        return np.sqrt( np.mean( (np.transpose(solution.y) - data) ** 2 ) )


    def fit(self):
        '''
        Estimate β and ɣ given integrated values from SIR model (cumulative sums).
        To fit the curve (thus getting β and ɣ) we must minimize the error using RMS.
        :return: Estimated beta and gamma
        '''
        optimal = minimize(
            self.loss_rms,              #Loss function
            [0.001, 0.001],             #Initial guess
            args = (self.data),
            method = 'L-BFGS-B',
            bounds = [(0.0001, 1.0), (0.0001, 1.0)]     #β, ɣ bounds
        )

        self.beta, self.gamma = optimal.x

        return self.beta, self.gamma




def main():

    italian_population = 500000

    import data_fetcher as df
    raw_data = df.fetch_data()
    data = df.process_data(raw_data, population=italian_population)
    inter = SIRInterpolation(data)
    beta, gamma = inter.fit()
    print(beta, gamma)


    print("Come si evolverà la situazione?")
    model = SIRPredict(italian_population, beta, gamma, len(data), data[-1])

    predicted_data = np.ndarray(shape=(100, 3))
    for i, _ in enumerate(predicted_data):
        predicted_data[i] = model.solve()

    plt.plot(data)
    plt.show()

    plt.plot(np.vstack((data, predicted_data)))
    plt.show()


if __name__ == "__main__":
    main()