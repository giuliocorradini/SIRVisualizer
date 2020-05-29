import cellpylib as cpl
import random
import numpy as np
import matplotlib.pyplot as plt

class VirusSpreadAutomata:

    SUSCEPTIBLE = 0
    INFECTIOUS = 1
    RECOVERED = 2

    def __init__(self, nb_size: int = 10, tr: float = 1/20, rr: float = 1/2):
        self.sandbox = cpl.init_simple2d(nb_size, nb_size)

        center = int(nb_size/2)
        self.sandbox[:, center, center] = self.INFECTIOUS

        self.population = nb_size * nb_size
        self.susceptible = self.population - 1
        self.infectious = 1
        self.recovered = 0

        self.transmission_rate = tr
        self.recover_rate = rr

        self.time = 1

    def virus_spread_rule(self, neighbourhood, c, t):
        cell = neighbourhood[1][1]

        if cell == self.SUSCEPTIBLE:
            infectious_neighs = np.count_nonzero(neighbourhood == self.INFECTIOUS)
            comp_risk = infectious_neighs * random.random()
            if comp_risk >= self.transmission_rate:
                return self.INFECTIOUS
            else:
                return self.SUSCEPTIBLE

        elif cell == self.INFECTIOUS:
            if random.random() >= self.recover_rate:
                return self.RECOVERED
            else:
                return self.INFECTIOUS

        else:
            return self.RECOVERED

    def evolve(self, ts: int = 1):
        self.sandbox = cpl.evolve2d(self.sandbox, timesteps=ts, neighbourhood='Moore',
                                    apply_rule=lambda n,c,t: self.virus_spread_rule(n,c,t))

    def get_sir(self):
        for nb_t in self.sandbox:
            susceptible = np.count_nonzero(nb_t == self.SUSCEPTIBLE)
            infectious = np.count_nonzero(nb_t == self.INFECTIOUS)
            recovered = self.population - susceptible - infectious

            yield (susceptible, infectious, recovered)


def main():
    vsa = VirusSpreadAutomata(100, 0.5573631610109397, 0.45333239503416245)
    vsa.evolve(90)

    cpl.plot2d_animate(vsa.sandbox)

    plt.plot(list(vsa.get_sir()))
    plt.legend(('Susceptible', 'Infectious', 'Recovered'))
    plt.show()


if __name__ == '__main__':
    main()