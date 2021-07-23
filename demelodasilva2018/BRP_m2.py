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


class BRP_m2:
    def __init__(self, bay, bugs_fixed=True, restricted=False, distinct=False):
        if restricted or distinct:
            assert not restricted or distinct
            bay.validate_distinct()

        self.model = model = Model()

        self.S = S = bay.n_stacks
        self.H = H = bay.n_tiers
        self.N = N = bay.n_blocks
        self.G = G = bay.p_max
        self.T = T = bay.brp_min_max()

        self.x = x = model.binary_var_dict(product(irange(0, T), irange(1, G), irange(1, S), irange(1, H)), name='x')
        self.y = y = model.binary_var_dict(product(irange(1, T), irange(1, G), irange(1, S), irange(1, H)), name='y')
        self.z = z = model.binary_var_dict(product(irange(1, T), irange(1, G), irange(1, S), irange(1, H)), name='z')
        self.k = k = model.binary_var_dict(product(irange(0, T), irange(1, G), irange(1, S), irange(1, H)), name='k')

        C = {(g, s, h): int(g == bay.pri[s - 1][h - 1]) for g, s, h in product(irange(1, G), irange(1, S), irange(1, H))}

        # objective
        model.minimize(model.sum(y[t, g, s, h] for t, g, s, h in product(irange(1, T), irange(1, G), irange(1, S), irange(1, H))))
        # (2)
        if bugs_fixed:
            model.add_constraints(x[0, g, s, h] + k[0, g, s, h] == C[g, s, h] for g, s, h in product(irange(1, G), irange(1, S), irange(1, H)))
        else:
            model.add_constraints(x[0, g, s, h] == C[g, s, h] for g, s, h in product(irange(1, G), irange(1, S), irange(1, H)))
            model.add_constraints(k[0, g, s, h] == 0 for g, s, h in product(irange(1, G), irange(1, S), irange(1, H)))
        # (4)
        model.add_constraints(model.sum(x[t, g, s, h] for g in irange(1, G)) >= model.sum(x[t, g, s, h + 1] for g in irange(1, G)) for t, s, h in product(irange(1, T), irange(1, S), irange(1, H - 1)))
        # (8)
        model.add_constraints(model.sum(y[t, g, s, h] for g, s, h in product(irange(1, G), irange(1, S), irange(1, H))) <= 1 for t in irange(1, T))
        # (9)
        model.add_constraints(model.sum(z[t, g, s, h] for g, s, h in product(irange(1, G), irange(1, S), irange(1, H))) <= 1 for t in irange(1, T))
        if bugs_fixed:
            # relocations are forced to be carried out as early as possible
            model.add_constraints(model.sum(z[t - 1, g, s, h] for g, s, h in product(irange(1, G), irange(1, S), irange(1, H))) >= model.sum(z[t, g, s, h] for g, s, h in product(irange(1, G), irange(1, S), irange(1, H))) for t in irange(2, T))
            # retrievals only occur after a relocation (except for time 0)
            model.add_constraints(model.sum(k[t, g, s, h] for g, s, h in product(irange(1, G), irange(1, S), irange(1, H))) <= N * model.sum(z[t, g, s, h] for g, s, h in product(irange(1, G), irange(1, S), irange(1, H))) for t in irange(1, T))
        # (15)
        model.add_constraints(model.sum(y[t, g, s, 1] for g in irange(1, G)) <= 1 - model.sum(x[t - 1, g, s, 1] for g in irange(1, G)) for t, s in product(irange(1, T), irange(1, S)))
        # (16)
        model.add_constraints(model.sum(y[t, g, s, h + 1] for g in irange(1, G)) <= model.sum(x[t - 1, g, s, h] - x[t - 1, g, s, h + 1] for g in irange(1, G)) for t, s, h in product(irange(1, T), irange(1, S), irange(1, H - 1)))
        # (17)
        model.add_constraints(model.sum(z[t, g, s, h] for g in irange(1, G)) <= model.sum(x[t - 1, g, s, h] - x[t - 1, g, s, h + 1] for g in irange(1, G)) for t, s, h in product(irange(1, T), irange(1, S), irange(1, H - 1)))
        # (20)
        model.add_constraints(model.sum(z[t, g, s, h] + y[t, g, s, h] for g, h in product(irange(1, G), irange(1, H))) <= 1 for t, s in product(irange(1, T), irange(1, S)))
        # (41)
        model.add_constraint(model.sum(x[T, g, s, h] for g, s, h in product(irange(1, G), irange(1, S), irange(1, H))) == 0)
        # (42)
        model.add_constraints(model.sum(x[t, g, s, h] + k[t, g, s, h] for g in irange(1, G)) <= 1 for t, s, h in product(irange(1, T), irange(1, S), irange(1, H)))
        # (43)
        model.add_constraints(model.sum(x[t - 1, g, s, h] for s, h in product(irange(1, S), irange(1, H))) == model.sum(x[t, g, s, h] + k[t, g, s, h] for s, h in product(irange(1, S), irange(1, H))) for t, g in product(irange(1, T), irange(1, G)))
        # (44)
        model.add_constraints(z[t, g, s, h] + k[t, g, s, h] == y[t, g, s, h] + x[t - 1, g, s, h] - x[t, g, s, h] for t, g, s, h in product(irange(1, T), irange(1, G), irange(1, S), irange(1, H)))

        if distinct:
            # (46)
            if bugs_fixed:
                model.add_constraints(model.sum(k[t, g, s, h] for s, h in product(irange(1, S), irange(1, H))) <= model.sum(k[u, g - 1, s, h] for u, s, h in product(irange(0, t), irange(1, S), irange(1, H))) for t, g in product(irange(0, T), irange(2, G)))
            else:
                model.add_constraints(model.sum(k[t, g, s, h] for s, h in product(irange(1, S), irange(1, H))) <= model.sum(k[u, g - 1, s, h] for u, s, h in product(irange(1, t), irange(1, S), irange(1, H))) for t, g in product(irange(1, T), irange(2, G)))
            # (47)
            if bugs_fixed:
                model.add_constraints(k[t, g, s, h] + model.sum(k[t, l, s, h + 1] for l in irange(g + 1, G)) <= 1 for t, g, s, h in product(irange(0, T), irange(1, G - 1), irange(1, S), irange(1, H - 1)))
            else:
                model.add_constraints(k[t, g, s, h] + model.sum(k[t, l, s, h + 1] for l in irange(g + 1, G)) <= 1 for t, g, s, h in product(irange(1, T), irange(1, G - 1), irange(1, S), irange(1, H - 1)))
            # (48)
            if bugs_fixed:
                model.add_constraints(k[t, g, s, h] <= 1 - model.sum(x[t, l, s, h - 1] for l in irange(1, g - 1)) for t, g, s, h in product(irange(0, T), irange(2, G), irange(1, S), irange(2, H)))
            else:
                model.add_constraints(k[t, g, s, h] <= 1 - model.sum(x[t, l, s, h - 1] for l in irange(1, g - 1)) for t, g, s, h in product(irange(1, T), irange(2, G), irange(1, S), irange(2, H)))
        else:
            # (49)
            if bugs_fixed:
                model.add_constraints(k[t, g, s, h] <= 1 - model.sum(x[t, l, s, i] + k[t, l, s, i] for l in irange(1, g - 1)) for t, g, s, h in product(irange(0, T), irange(2, G), irange(1, S), irange(2, H)) for i in irange(1, h - 1))
            else:
                model.add_constraints(k[t, g, s, h] <= 1 - model.sum(x[t, l, s, i] for l in irange(1, g - 1)) for t, g, s, h in product(irange(1, T), irange(2, G), irange(1, S), irange(2, H)) for i in irange(1, h - 1))
            # (50)
            if bugs_fixed:
                model.add_constraints(k[t, g, s, h] <= 1 - model.sum(x[t, l, p, i] for l in irange(1, g - 1)) for t, g, s, h in product(irange(0, T), irange(2, G), irange(1, S), irange(1, H)) for p in irange(1, S) if p != s for i in irange(1, H))
            else:
                model.add_constraints(k[t, g, s, h] <= 1 - model.sum(x[t, l, p, i] for l in irange(1, g - 1)) for t, g, s, h in product(irange(1, T), irange(2, G), irange(1, S), irange(1, H)) for p in irange(1, S) if p != s for i in irange(1, H))

        if restricted:
            # (52)
            if bugs_fixed:
                model.add_constraints(model.sum(z[t, j, s, h] for j, h in product(irange(g + 1, G), irange(1, H))) <= model.sum(k[u, g, r, l] for u, r, l in product(irange(1, t - 1), irange(1, S), irange(1, H))) + model.sum(x[t - 1, i, s, h] for i, h in product(irange(1, g), irange(1, H))) for t, g, s in product(irange(1, T), irange(1, G - 1), irange(1, S)))
            else:
                model.add_constraints(model.sum(z[t, j, s, h] for j, h in product(irange(g + 1, G), irange(1, H))) <= model.sum(k[u, g, r, l] for u, r, l in product(irange(1, t), irange(1, S), irange(1, H))) + model.sum(x[t, i, s, h] for i, h in product(irange(1, g), irange(1, H))) for t, g, s in product(irange(1, T), irange(1, G - 1), irange(1, S)))

    def get_bays(self):
        before = {}
        after = {}
        for t in irange(0, self.T):
            conf_before = [[None] * self.H for _ in range(self.S)]
            conf_after = [[None] * self.H for _ in range(self.S)]
            for g, s, h in product(irange(1, self.G), irange(1, self.S), irange(1, self.H)):
                if round(self.x[t, g, s, h].solution_value) + round(self.k[t, g, s, h].solution_value) == 1:
                    conf_before[s - 1][h - 1] = g
                if round(self.x[t, g, s, h].solution_value) == 1:
                    conf_after[s - 1][h - 1] = g
            before[t] = Bay(self.S, self.H, conf_before)
            after[t] = Bay(self.S, self.H, conf_after)
        return before, after

    def get_n_relos(self):
        return self.model.objective_value


def test():
    conf = [[4, 1], [2], [1, 3, 4]]
    bay = Bay(3, 3, conf)
    print('bay')
    print(bay)

    brp_m2 = BRP_m2(bay)
    if brp_m2.model.solve():
        print()
        print('n_relos = {}'.format(brp_m2.get_n_relos()))
        before, after = brp_m2.get_bays()
        for t in irange(0, brp_m2.T):
            print('t = {} (before retrievals)'.format(t))
            print(before[t])
            print('t = {} (after retrievals)'.format(t))
            print(after[t])


if __name__ == '__main__':
    test()
