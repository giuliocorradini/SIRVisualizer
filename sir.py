import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize


class SIRPredict:
    '''
    Represents the SIR model facilities to compute and predict the spread of
    a virus given the value of β and ɣ
    '''

    def __init__(self, population: int, beta: float, gamma: float):
        self.beta = beta
        self.gamma = gamma

        self.n = population

        self.reset_model()

    def reset_model(self):
        '''
        Initialize and resets model for computation.
        Susceptible = population - 1
        Infected = 1 (index case)
        '''
        self.s = self.n - 1
        self.i = 1
        self.r = 0

    def solve(self):
        s_next = -(self.beta * self.s * self.i / self.n)
        i_next = (self.beta * self.s * self.i / self.n) - self.gamma * self.i
        r_next = self.gamma * self.i

        self.s = s_next
        self.i = i_next
        self.r = r_next

        return self.s, self.i, self.r

    def predict(self, day: int):
        '''
        Predict epidemic spread given the day, in the discrete domain.
        :param day: Day from the beginning of spread
        :return: Tuple with SIR values
        '''
        for _ in range(day):
            self.solve()

        return self.s, self.i, self.r

class SIRInterpolation:
    def __init__(self, sir_values: np.ndarray):
        '''
        Initializes the SIR model given an array of SIR functions first discrete derivative.
        :param sir_values: Array with SIR values.
        '''
        self.data = sir_values
        self.n = sir_values[0][0] + 1

    @classmethod
    def loss_rms(cls, point, data):
        '''
        Tries given
        :param point:
        :param data:
        :return: float
        '''

        size = len(data)
        beta, gamma = point

        def next_dt_SIR(time, sir_data):
            s, i, r = sir_data
            return [-beta*s*i, beta*s*i-gamma*i, gamma*i]

        solution = solve_ivp(next_dt_SIR,
                             (0, size),
                             [S_0, I_0, R_0],
                             t_eval=np.arange(0, size, 1),
                             vectorized=True
                             )
        return np.sqrt(np.mean((solution.y[1]-data)**2))


    def interpolate(self):
        '''
        Estimate β and ɣ given integrated values from SIR model (cumulative sums).
        To fit the curve (thus getting β and ɣ) we must minimize the error using RMS.
        :return:
        '''
        optimal = minimize(
            SIRInterpolation.loss_rms,  #Loss function
            [0.001, 0.001],             #Initial guess
            args = (self.data),
            method = 'L-BFGS-B',
            bounds = [(0.00000001, 0.4), (0.00000001, 0.4)]     #β, ɣ bounds
        )

        self.beta, self.gamma = optimal.x

        return self.beta, self.gamma




def main():
    model = SIRPredict(1000, 0.2, 0.1)
    plt.plot([model.solve() for _ in range(200)])
    plt.show()


if __name__ == "__main__":
    main()