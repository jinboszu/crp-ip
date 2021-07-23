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


class CRP_I:
    def __init__(self, bay):
        bay.validate_full_distinct()
        self.model = model = Model()

        self.S = S = bay.n_stacks
        self.T = T = bay.n_tiers
        self.C = C = bay.n_blocks
        self.N = N = C - S + 1

        self.a = a = model.binary_var_dict(((n, c, d) for n in irange(1, N) for c, d in product(irange(n, C + S), irange(n, C))), name='a')
        self.b = b = model.binary_var_dict(irange(N + 1, C), name='b')

        stack = {bay.pri[s][t]: s + 1 for s in range(bay.n_stacks) for t in range(bay.h[s])}
        tier = {bay.pri[s][t]: t + 1 for s in range(bay.n_stacks) for t in range(bay.h[s])}
        A1 = {(c, d): int(stack[c] == stack[d] and tier[c] < tier[d]) for c, d in product(irange(1, C), irange(1, C))}
        A1.update({(C + s, d): int(s == stack[d]) for s, d in product(irange(1, S), irange(1, C))})

        # objective
        model.minimize(model.sum(a[n, n, d] for n in irange(1, N - 1) for d in irange(n + 1, C)) + model.sum(b[d] for d in irange(N + 1, C)))
        # (1)
        model.add_constraints(a[1, c, d] == A1[c, d] for c, d in product(irange(1, C + S), irange(1, C)))
        # (2)
        model.add_constraints(model.sum(a[n, C + s, c] for s in irange(1, S)) == 1 for n in irange(2, N) for c in irange(n, C))
        # (3)
        model.add_constraints(a[n, c, c] == 0 for n in irange(2, N) for c in irange(n, C))
        # (4)
        model.add_constraints(a[n, c, d] + a[n, d, c] <= 1 for n in irange(2, N) for c, d in product(irange(n, C), irange(n, C)) if c != d)
        # (5)
        model.add_constraints(a[n, c, d] + a[n, d, c] >= a[n, C + s, c] + a[n, C + s, d] - 1 for n, s in product(irange(2, N), irange(1, S)) for c, d in product(irange(n, C), irange(n, C)) if c != d)
        # (6)
        model.add_constraints(a[n, c, d] + a[n, d, c] <= 2 - a[n, C + s, c] - model.sum(a[n, C + r, d] for r in irange(1, S) if r != s) for n, s in product(irange(2, N), irange(1, S)) for c, d in product(irange(n, C), irange(n, C)) if c != d)
        # (7)
        model.add(model.sum(a[n, C + s, d] for d in irange(n, C)) <= T for n, s in product(irange(2, N), irange(1, S)))
        # (8)
        model.add_constraints(b[d] >= a[N, c, d] for d in irange(N + 1, C) for c in irange(N, d - 1))
        # (9)
        model.add_constraints(a[n + 1, d, c] <= a[n, d, c] + a[n, n, c] for n in irange(1, N - 1) for c, d in product(irange(n + 1, C), irange(n + 1, C + S)) if c != d)
        # (10)
        model.add_constraints(a[n + 1, d, c] >= a[n, d, c] - a[n, n, c] for n in irange(1, N - 1) for c, d in product(irange(n + 1, C), irange(n + 1, C + S)) if c != d)
        # (11)
        model.add_constraints(a[n, n, c] + a[n, C + s, c] + a[n + 1, C + s, c] <= 2 for n in irange(1, N - 1) for c, s in product(irange(n + 1, C), irange(1, S)))
        # (12)
        model.add_constraints(a[n, n, c] + a[n, n, d] + a[n, c, d] + a[n + 1, c, d] <= 3 for n in irange(1, N - 1) for c, d in product(irange(n + 1, C), irange(n + 1, C)) if c != d)

    def get_bays(self):
        bays = {}
        for n in irange(1, self.N):
            conf = [[None] * self.T for _ in range(self.S)]
            for c in irange(n, self.C):
                stack = next(s - 1 for s in irange(1, self.S) if round(self.a[n, self.C + s, c].solution_value) == 1)
                tier = sum(round(self.a[n, d, c].solution_value) for d in irange(n, self.C))
                conf[stack][tier] = c
            bays[n] = Bay(self.S, self.T, conf)
        return bays

    def get_n_relos(self):
        return self.model.objective_value


def test():
    conf = [[4], [3, 1], [2, 5, 6]]
    bay = Bay(3, 3, conf)
    print('bay')
    print(bay)

    crp_i = CRP_I(bay)
    if crp_i.model.solve():
        print()
        print('n_relos = {}'.format(crp_i.get_n_relos()))
        bays = crp_i.get_bays()
        for n in irange(1, crp_i.N):
            print('n = {}'.format(n))
            print(bays[n])


if __name__ == '__main__':
    test()
