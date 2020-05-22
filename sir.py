import matplotlib.pyplot as plt


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

    def ds(self):
        self.s += -(self.beta * self.s * self.i / self.n)
        return self.s

    def di(self):
        self.i += (self.beta * self.s * self.i / self.n) - self.gamma * self.i
        return self.i

    def dr(self):
        self.r += self.gamma * self.i
        return self.r

    def solve(self):
        self.ds()
        self.di()
        self.dr()

        return self.s, self.i, self.r

    def predict(self, day: int):
        '''
        Predict epidemic spread given the day, in the discrete domain.
        :param day: Day from the beginning of spread
        :return: Tuple with SIR values
        '''
        for _ in range(day):
            self.ds()
            self.di()
            self.dr()

        return self.s, self.i, self.r


def main():
    model = SIR(1000, 0.2, 0.1)
    plt.plot([model.solve() for _ in range(200)])
    plt.show()


if __name__ == "__main__":
    main()