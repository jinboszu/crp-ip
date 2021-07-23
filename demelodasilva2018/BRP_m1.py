# Copyright (c) 2021 Bo Jin <jinbostar@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


from itertools import product

from docplex.mp.model import Model

from bay import Bay
from common import irange


class BRP_m1:
    def __init__(self, bay, restricted_distinct=False):
        if restricted_distinct:
            bay.validate_distinct()
        self.model = model = Model()

        self.S = S = bay.n_stacks
        self.H = H = bay.n_tiers
        self.G = G = bay.p_max
        self.N = N = bay.n_blocks
        self.T = T = N + bay.brp_min_max()

        self.x = x = model.binary_var_dict(product(irange(0, T), irange(1, G), irange(1, S), irange(1, H)), name='x')
        self.y = y = model.binary_var_dict(product(irange(1, T), irange(1, G), irange(1, S), irange(1, H)), name='y')
        self.z = z = model.binary_var_dict(product(irange(1, T), irange(1, G), irange(1, S), irange(1, H)), name='z')
        self.k = k = model.binary_var_dict(product(irange(1, T), irange(1, G), irange(1, N)), name='k')
        self.w = w = model.binary_var_dict(product(irange(0, T), irange(1, G), irange(1, N)), name='w')

        C = {(g, s, h): int(g == bay.pri[s - 1][h - 1]) for g, s, h in product(irange(1, G), irange(1, S), irange(1, H))}
        values = sorted(bay.values)
        Q = {(g, n): int(g == values[n - 1]) for g, n in product(irange(1, G), irange(1, N))}

        # objective
        model.minimize(model.sum(y[t, g, s, h] for t, g, s, h in product(irange(1, T), irange(1, G), irange(1, S), irange(1, H))))
        # (2)
        model.add_constraints(x[0, g, s, h] == C[g, s, h] for g, s, h in product(irange(1, G), irange(1, S), irange(1, H)))
        # (4)
        model.add_constraints(model.sum(x[t, g, s, h] for g in irange(1, G)) >= model.sum(x[t, g, s, h + 1] for g in irange(1, G)) for t, s, h in product(irange(1, T), irange(1, S), irange(1, H - 1)))
        # (8)
        model.add_constraints(model.sum(y[t, g, s, h] for g, s, h in product(irange(1, G), irange(1, S), irange(1, H))) <= 1 for t in irange(1, T))
        # (9)
        model.add_constraints(model.sum(z[t, g, s, h] for g, s, h in product(irange(1, G), irange(1, S), irange(1, H))) <= 1 for t in irange(1, T))
        # (11)
        model.add_constraints(model.sum(z[t - 1, g, s, h] for g, s, h in product(irange(1, G), irange(1, S), irange(1, H))) >= model.sum(z[t, g, s, h] for g, s, h in product(irange(1, G), irange(1, S), irange(1, H))) for t in irange(2, T))
        # (15)
        model.add_constraints(model.sum(y[t, g, s, 1] for g in irange(1, G)) <= 1 - model.sum(x[t - 1, g, s, 1] for g in irange(1, G)) for t, s in product(irange(1, T), irange(1, S)))
        # (16)
        model.add_constraints(model.sum(y[t, g, s, h + 1] for g in irange(1, G)) <= model.sum(x[t - 1, g, s, h] - x[t - 1, g, s, h + 1] for g in irange(1, G)) for t, s, h in product(irange(1, T), irange(1, S), irange(1, H - 1)))
        # (17)
        model.add_constraints(model.sum(z[t, g, s, h] for g in irange(1, G)) <= model.sum(x[t - 1, g, s, h] - x[t - 1, g, s, h + 1] for g in irange(1, G)) for t, s, h in product(irange(1, T), irange(1, S), irange(1, H - 1)))
        # (19)
        model.add_constraints(model.sum(x[t, g, s, h] + z[t, g, s, h] for g in irange(1, G)) <= 1 for t, s, h in product(irange(1, T), irange(1, S), irange(1, H)))
        # (25)
        model.add_constraint(model.sum(x[T, g, s, h] for g, s, h in product(irange(1, G), irange(1, S), irange(1, H))) == 0)
        # (26)
        model.add_constraints(w[0, g, n] == 0 for g, n in product(irange(1, G), irange(1, N)))
        # (27)
        model.add_constraints(w[T, g, n] == Q[g, n] for g, n in product(irange(1, G), irange(1, N)))
        # (28)
        model.add_constraints(model.sum(w[t, g, n] for g in irange(1, G)) <= 1 for t, n in product(irange(1, T), irange(1, N)))
        # (29)
        model.add_constraints(model.sum(k[t, g, n] for g, n in product(irange(1, G), irange(1, N))) <= 1 for t in irange(1, T))
        # (30)
        model.add_constraints(model.sum(w[t, g, n] for g in irange(1, G)) >= model.sum(w[t, g, n + 1] for g in irange(1, G)) for t, n in product(irange(1, T), irange(1, N - 1)))
        # (31)
        model.add_constraints(model.sum(w[t, g, n] for n in irange(1, N)) == model.sum(w[t - 1, g, n] + k[t, g, n] for n in irange(1, N)) for t, g in product(irange(1, T), irange(1, G)))
        # (32)
        model.add_constraints(w[t, g, n] == w[t - 1, g, n] + k[t, g, n] for t, g, n in product(irange(1, T), irange(1, G), irange(1, N)))
        # (33)
        model.add_constraints(model.sum(x[t - 1, g, s, h] for s, h in product(irange(1, S), irange(1, H))) == model.sum(k[t, g, n] for n in irange(1, N)) + model.sum(x[t, g, s, h] for s, h in product(irange(1, S), irange(1, H))) for t, g in product(irange(1, T), irange(1, G)))
        # (34)
        model.add_constraints(x[t, g, s, h] + z[t, g, s, h] >= x[t - 1, g, s, h] + y[t, g, s, h] for t, g, s, h in product(irange(1, T), irange(1, G), irange(1, S), irange(1, H)))
        # (35)
        model.add_constraints(model.sum(y[t, g, s, h] for s, h in product(irange(1, S), irange(1, H))) + model.sum(k[t, g, n] for n in irange(1, N)) == model.sum(z[t, g, s, h] for s, h in product(irange(1, S), irange(1, H))) for t, g in product(irange(1, T), irange(1, G)))
        # (36)
        model.add_constraints(model.sum(z[t, g, s, h] for g, s, h in product(irange(1, G), irange(1, S), irange(1, H))) == 1 for t in irange(1, N))
        if restricted_distinct:
            # (51)
            model.add_constraints(model.sum(z[t, j, s, h] for j, h in product(irange(g + 1, G), irange(1, H))) <= w[t - 1, g, g] + model.sum(x[t, i, s, h] for i, h in product(irange(1, g), irange(1, H))) for t, g, s in product(irange(1, T), irange(1, G - 1), irange(1, S)))

    def get_bays(self):
        bays = {}
        for t in irange(0, self.T):
            conf = [[None] * self.H for _ in range(self.S)]
            for g, s, h in product(irange(1, self.G), irange(1, self.S), irange(1, self.H)):
                if round(self.x[t, g, s, h].solution_value) == 1:
                    conf[s - 1][h - 1] = g
            bays[t] = Bay(self.S, self.H, conf)
        return bays

    def get_n_relos(self):
        return self.model.objective_value


def test1():
    conf = [[2], [2], [1, 3, 4]]
    bay = Bay(3, 3, conf)
    print('bay')
    print(bay)

    brp_m1 = BRP_m1(bay)
    if brp_m1.model.solve():
        print()
        print('n_relos = {}'.format(brp_m1.get_n_relos()))
        bays = brp_m1.get_bays()
        for t in irange(0, brp_m1.T):
            print('t = {}'.format(t))
            print(bays[t])


def test2():
    conf = [[3], [2], [1, 4, 5]]
    bay = Bay(3, 3, conf)
    print('bay')
    print(bay)

    brp_m1 = BRP_m1(bay, restricted_distinct=True)
    if brp_m1.model.solve():
        print()
        print('n_relos = {}'.format(brp_m1.get_n_relos()))
        bays = brp_m1.get_bays()
        for t in irange(0, brp_m1.T):
            print('t = {}'.format(t))
            print(bays[t])


if __name__ == '__main__':
    test1()
    print()
    print()
    test2()
