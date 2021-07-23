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


class BRP_III:
    def __init__(self, bay):
        bay.validate_full_distinct()
        self.model = model = Model()

        self.S = S = bay.n_stacks
        self.mxHeight = mxHeight = bay.n_tiers
        self.C = C = bay.n_blocks
        self.W = W = C + bay.brp_min_max()

        self.X3 = X3 = model.continuous_var_dict(product(irange(1, C), irange(1, S), irange(1, W + 1)), lb=0, ub=1, name='X')
        self.B2 = B2 = model.continuous_var_dict(product(irange(1, C), irange(1, W + 1)), lb=0, ub=mxHeight, name='B')
        self.M2 = M2 = model.binary_var_dict(product(irange(1, C), irange(1, W)), name='M')
        self.C2 = C2 = model.binary_var_dict(product(irange(1, C), irange(1, W)), name='C')
        self.F2 = F2 = model.binary_var_dict(product(irange(1, C), irange(1, W)), name='F')
        self.T2 = T2 = model.binary_var_dict(product(irange(1, C), irange(1, W)), name='T')
        self.R2 = R2 = model.binary_var_dict(product(irange(1, S), irange(1, W)), name='R')
        self.P2 = P2 = model.binary_var_dict(product(irange(1, S), irange(1, W)), name='P')
        self.R3 = R3 = model.continuous_var_dict(product(irange(1, C), irange(1, S), irange(1, W)), lb=0, ub=1, name='R')
        self.P3 = P3 = model.continuous_var_dict(product(irange(1, C), irange(1, S), irange(1, W)), lb=0, ub=1, name='P')

        initialStack = {bay.pri[s][t]: s + 1 for s in range(bay.n_stacks) for t in range(bay.h[s])}
        initialSetup = {(c, s): int(initialStack[c] == s) for c, s in product(irange(1, C), irange(1, S))}
        initialBury = {bay.pri[s][t]: bay.h[s] - t for s in range(bay.n_stacks) for t in range(bay.h[s])}

        # objective
        model.minimize(model.sum(t * T2[C, t] for t in irange(1, W)))
        # (1a)
        model.add_constraints(R3[c, s, t] <= R2[s, t] for c, s, t in product(irange(1, C), irange(1, S), irange(1, W)))
        # (1b)
        model.add_constraints(R3[c, s, t] <= M2[c, t] for c, s, t in product(irange(1, C), irange(1, S), irange(1, W)))
        # (1c)
        model.add_constraints(R3[c, s, t] >= R2[s, t] + M2[c, t] - 1 for c, s, t in product(irange(1, C), irange(1, S), irange(1, W)))
        # (1e)
        model.add_constraints(model.sum(R3[c, s, t] for c in irange(1, C)) == R2[s, t] for s, t in product(irange(1, S), irange(1, W)))
        # (1f)
        model.add_constraints(model.sum(R3[c, s, t] for s in irange(1, S)) == M2[c, t] for c, t in product(irange(1, C), irange(1, W)))
        # (2a)
        model.add_constraints(P3[c, s, t] <= P2[s, t] for c, s, t in product(irange(1, C), irange(1, S), irange(1, W)))
        # (2b)
        model.add_constraints(P3[c, s, t] <= M2[c, t] for c, s, t in product(irange(1, C), irange(1, S), irange(1, W)))
        # (2c)
        model.add_constraints(P3[c, s, t] >= P2[s, t] + M2[c, t] - 1 for c, s, t in product(irange(1, C), irange(1, S), irange(1, W)))
        # (2e)
        model.add_constraints(model.sum(P3[c, s, t] for c in irange(1, C)) == P2[s, t] for s, t in product(irange(1, S), irange(1, W)))
        # (2f)
        model.add_constraints(model.sum(P3[c, s, t] for s in irange(1, S)) <= M2[c, t] for c, t in product(irange(1, C), irange(1, W)))
        # (3a)
        model.add_constraints(X3[c, s, 1] == initialSetup[c, s] for c, s in product(irange(1, C), irange(1, S)))
        # (3b)
        model.add_constraints(X3[c, s, t + 1] == X3[c, s, t] + P3[c, s, t] - R3[c, s, t] for c, s, t in product(irange(1, C), irange(1, S), irange(1, W)))
        # (4a)
        model.add_constraints(B2[c, 1] == initialBury[c] for c in irange(1, C))
        # (4b)
        model.add_constraints(B2[c, t + 1] == B2[c, t] + F2[c, t] - C2[c, t] for c, t in product(irange(1, C), irange(1, W)))
        # (5)
        model.add_constraints(model.sum(X3[c, s, t] for s in irange(1, S)) <= 1 for c, t in product(irange(1, C), irange(1, W + 1)))
        # (6)
        model.add_constraints(model.sum(X3[c, s, t] for c in irange(1, C)) <= mxHeight for s, t in product(irange(1, S), irange(1, W + 1)))
        # (7*) safety constraints not included
        # model.add_constraints(model.sum(X3[c, s, t] for c in irange(1, C)) <= 2 + model.sum(X3[c, s + 1, t] for c in irange(1, C)) for s, t in product(irange(1, S - 1), irange(1, W + 1)))
        # (8*) safety constraints not included
        # model.add_constraints(model.sum(X3[c, s, t] for c in irange(1, C)) >= model.sum(X3[c, s + 1, t] for c in irange(1, C)) - 2 for s, t in product(irange(1, S - 1), irange(1, W + 1)))
        # (10)
        model.add_constraints(model.sum(M2[c, t] for t in irange(1, W)) >= 1 for c in irange(1, C))
        # (11)
        model.add_constraints(model.sum(M2[c, t] for c in irange(1, C)) <= 1 for t in irange(1, W))
        # (12)
        model.add_constraints(model.sum(M2[c, t] for c in irange(1, C)) == 1 for t in irange(1, C))
        # (13)
        model.add_constraints(model.sum(C2[c, t] for c in irange(1, C)) <= mxHeight for t in irange(1, W))
        # (14)
        model.add_constraints(model.sum(F2[c, t] for c in irange(1, C)) <= mxHeight for t in irange(1, W))
        # (15)
        model.add_constraints(model.sum(C2[c, t] - F2[c, t] for t in irange(1, W)) == initialBury[c] for c in irange(1, C))
        # (16)
        model.add_constraints(model.sum(T2[c, t] for c in irange(1, C)) <= 1 for t in irange(1, W))
        # (17)
        model.add_constraints(model.sum(T2[c, t] for t in irange(1, W)) == 1 for c in irange(1, C))
        # (18)
        model.add_constraints(model.sum(t * T2[c + 1, t] for t in irange(1, W)) >= 1 + model.sum(t * T2[c, t] for t in irange(1, W)) for c in irange(1, C - 1))
        # (19)
        model.add_constraints(model.sum(R2[s, t] for s in irange(1, S)) <= 1 for t in irange(1, W))
        # (20)
        model.add_constraints(model.sum(R2[s, t] for s in irange(1, S)) == 1 for t in irange(1, C))
        # (21)
        model.add_constraints(model.sum(P2[s, t] for s in irange(1, S)) <= model.sum(R2[s, t] for s in irange(1, S)) for t in irange(1, W))
        # (22)
        model.add_constraints(P2[s, t] + R2[s, t] <= 1 for s, t in product(irange(1, S), irange(1, W)))
        # (23)
        model.add_constraints(model.sum(P3[c, s, t] for s in irange(1, S)) <= model.sum(R3[c, s, t] for s in irange(1, S)) for c, t in product(irange(1, C), irange(1, W)))
        # (24)
        model.add_constraints(R3[c, s, t] <= X3[c, s, t] for c, s, t in product(irange(1, C), irange(1, S), irange(1, W)))
        # (25)
        model.add_constraints((1 - M2[c, t]) * (mxHeight - 1) >= B2[c, t] - 1 for c, t in product(irange(1, C), irange(1, W)))
        # (26)
        model.add_constraints(T2[c, t] <= M2[c, t] for c, t in product(irange(1, C), irange(1, W)))
        # (27)
        model.add_constraints(T2[c, t] == model.sum(R3[c, s, t] for s in irange(1, S)) - model.sum(P3[c, s, t] for s in irange(1, S)) for c, t in product(irange(1, C), irange(1, W)))
        # (28a)
        model.add_constraints(X3[c, s, t] >= C2[c, t] + R2[s, t] - 1 for c, s, t in product(irange(1, C), irange(1, S), irange(1, W)))
        # (28b)
        model.add_constraints(R2[s, t] >= X3[c, s, t] + C2[c, t] - 1 for c, s, t in product(irange(1, C), irange(1, S), irange(1, W)))
        # (28c)
        model.add_constraints(C2[c, t] >= R2[s, t] + X3[c, s, t] - 1 for c, s, t in product(irange(1, C), irange(1, S), irange(1, W)))
        # (29a)
        model.add_constraints(X3[c, s, t + 1] >= F2[c, t] + P2[s, t] - 1 for c, s, t in product(irange(1, C), irange(1, S), irange(1, W)))
        # (29b)
        model.add_constraints(P2[s, t] >= X3[c, s, t + 1] + F2[c, t] - 1 for c, s, t in product(irange(1, C), irange(1, S), irange(1, W)))
        # (29c)
        model.add_constraints(F2[c, t] >= P2[s, t] + X3[c, s, t + 1] - 1 for c, s, t in product(irange(1, C), irange(1, S), irange(1, W)))

    def get_bays(self):
        bays = {}
        for t in irange(1, self.W + 1):
            conf = [[None] * self.mxHeight for _ in range(self.S)]
            h = [0] * self.S
            for c, s in product(irange(1, self.C), irange(1, self.S)):
                if round(self.X3[c, s, t].solution_value) == 1:
                    conf[s - 1][round(self.B2[c, t].solution_value) - 1] = c
                    h[s - 1] += 1
            conf = [list(reversed(conf[s - 1][:h[s - 1]])) for s in irange(1, self.S)]
            bays[t] = Bay(self.S, self.mxHeight, conf)
        return bays

    def get_n_relos(self):
        return sum(round(self.P2[s, t].solution_value) for s, t in product(irange(1, self.S), irange(1, self.W)))


def test():
    conf = [[1, 3, 4], [5], [2]]
    bay = Bay(3, 3, conf)
    print('bay')
    print(bay)

    brp_iii = BRP_III(bay)
    if brp_iii.model.solve():
        print()
        print('n_relos = {}'.format(brp_iii.get_n_relos()))
        bays = brp_iii.get_bays()
        for t in irange(1, brp_iii.W + 1):
            print('t = {}'.format(t))
            print(bays[t])


if __name__ == '__main__':
    test()
