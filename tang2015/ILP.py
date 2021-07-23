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


class MRIP:
    def __init__(self, bay, bug_fixed=True):
        bay.validate_full_distinct()
        self.model = model = Model()

        self.C = C = bay.n_stacks
        self.P = P = bay.n_tiers
        self.S = S = bay.n_blocks

        self.x = x = model.binary_var_dict(((s, i, c, p) for s in irange(1, S) for i in irange(s, S) for c in irange(1, C) for p in irange(1, P)), name='x')
        self.y = y = model.binary_var_dict(((s, i) for s in irange(1, S - 1) for i in irange(s + 1, S)), name='y')
        self.w = w = model.binary_var_dict(((s, i, j) for s in irange(1, S - 1) for i in irange(s + 1, S) for j in irange(s + 1, S) if i != j), name='w')

        X1 = {(i, c, p): int(i == bay.pri[c - 1][p - 1]) for i, c, p in product(irange(1, S), irange(1, C), irange(1, P))}
        Q = {bay.pri[s][t]: bay.qlt[s][t] for s in range(bay.n_stacks) for t in range(bay.h[s])}

        # objective
        model.minimize(model.sum(y[s, i] for s in irange(1, S - 1) for i in irange(s + 1, S)))
        # (2)
        model.add_constraints((1 - model.sum(x[s, s, c, p] for p in irange(1, P))) * P + y[s, i] >= (model.sum(p * x[s, i, c, p] for p in irange(1, P)) - model.sum(p * x[s, s, c, p] for p in irange(1, P))) / P for s in irange(1, S - 1) for i in irange(s + 1, S) for c in irange(1, C))
        # (3)
        model.add_constraints((model.sum(p * x[s, s, c, p] for p in irange(1, P)) - model.sum(p * x[s, i, c, p] for p in irange(1, P))) / P <= 1 - y[s, i] for s in irange(1, S - 1) for i in irange(s + 1, S) for c in irange(1, C))
        # (4)
        model.add_constraints(model.sum(x[s, i, c, p] for c in irange(1, C) for p in irange(1, P)) == 1 for s in irange(1, S) for i in irange(s, S))
        # (5)
        model.add_constraints(model.sum(x[s, i, c, p] for i in irange(s, S)) <= 1 for s in irange(1, S) for c in irange(1, C) for p in irange(1, P))
        # (6)
        model.add_constraints(model.sum(x[s, i, c, p] for i in irange(s, S)) <= model.sum(x[s, i, c, p - 1] for i in irange(s, S)) for s in irange(1, S) for c in irange(1, C) for p in irange(2, P))
        # (7)
        model.add_constraints(model.sum(x[s + 1, i, c, p] for p in irange(1, P)) <= 2 - y[s, i] - model.sum(x[s, s, c, p] for p in irange(1, P)) for s in irange(1, S - 1) for i in irange(s + 1, S) for c in irange(1, C))
        # (8)
        model.add_constraints(2 - y[s, i] - y[s, j] + w[s, i, j] >= (model.sum(p * x[s, j, c, p] for c in irange(1, C) for p in irange(1, P)) - model.sum(p * x[s, i, c, p] for c in irange(1, C) for p in irange(1, P))) / P for s in irange(1, S - 1) for i in irange(s + 1, S) for j in irange(s + 1, S) if i != j)
        # (9)
        model.add_constraints(y[s, i] + y[s, j] + w[s, i, j] <= 3 + (model.sum(p * x[s, j, c, p] for c in irange(1, C) for p in irange(1, P)) - model.sum(p * x[s, i, c, p] for c in irange(1, C) for p in irange(1, P))) / P for s in irange(1, S - 1) for i in irange(s + 1, S) for j in irange(s + 1, S) if i != j)
        # (10)
        model.add_constraints(w[s, i, j] <= y[s, i] for s in irange(1, S - 1) for i in irange(s + 1, S) for j in irange(s + 1, S) if i != j)
        # (11)
        model.add_constraints(w[s, i, j] <= y[s, j] for s in irange(1, S - 1) for i in irange(s + 1, S) for j in irange(s + 1, S) if i != j)
        # (12)
        model.add_constraints(model.sum(p * x[s + 1, i, c, p] for p in irange(1, P)) - model.sum(p * x[s + 1, j, c, p] for p in irange(1, P)) >= - P * (1 - w[s, i, j]) - P * (1 - y[s, i]) - P * (1 - y[s, j]) - P * (1 - model.sum(x[s + 1, i, c, p] for p in irange(1, P))) for s in irange(1, S - 1) for i in irange(s + 1, S) for j in irange(s + 1, S) if i != j for c in irange(1, C))
        # (13)
        model.add_constraints(x[s + 1, i, c, p] - x[s, i, c, p] >= - y[s, i] for s in irange(1, S - 1) for i in irange(s + 1, S) for c in irange(1, C) for p in irange(1, P))
        # (25)
        model.add_constraints(x[s, i, c, p] - x[s + 1, i, c, p] >= - y[s, i] for s in irange(1, S - 1) for i in irange(s + 1, S) for c in irange(1, C) for p in irange(1, P))
        # (15)
        if bug_fixed:
            model.add_constraints(x[1, i, c, p] == X1[i, c, p] for i in irange(1, S) for c in irange(1, C) for p in irange(1, P))
        else:
            model.add_constraints(x[1, i, c, p] == X1[i, c, p] for i in irange(2, S) for c in irange(1, C) for p in irange(1, P))
        # (27)
        model.add_constraints(x[s, i, c, p] == X1[i, c, p] for i in irange(2, S) for s in irange(2, Q[i]) for c in irange(1, C) for p in irange(1, P))

    def get_bays(self):
        bays = {}
        for s in irange(1, self.S):
            conf = [[None] * self.P for _ in range(self.C)]
            for i, c, p in product(irange(s, self.S), irange(1, self.C), irange(1, self.P)):
                if round(self.x[s, i, c, p].solution_value) == 1:
                    conf[c - 1][p - 1] = i
            bays[s] = Bay(self.C, self.P, conf)
        return bays

    def get_n_relos(self):
        return self.model.objective_value


def test():
    conf = [[5, 3, 6], [4, 2], [7, 1]]
    bay = Bay(3, 3, conf)
    print('bay')
    print(bay)

    mrip = MRIP(bay)
    if mrip.model.solve():
        print()
        print('n_relos = {}'.format(mrip.get_n_relos()))
        bays = mrip.get_bays()
        for s in irange(1, mrip.S):
            print('s = {}'.format(s))
            print(bays[s])

    mrip = MRIP(bay, False)
    if mrip.model.solve():
        print()
        print('# bug unfixed #')
        print('n_relos = {}'.format(mrip.get_n_relos()))
        bays = mrip.get_bays()
        for s in irange(1, mrip.S):
            print('s = {}'.format(s))
            print(bays[s])


if __name__ == '__main__':
    test()
