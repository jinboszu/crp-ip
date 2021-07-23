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


class BRP_II_C:
    def __init__(self, bay):
        bay.validate_full_distinct()
        self.model = model = Model()

        self.W = W = bay.n_stacks
        self.H = H = bay.n_tiers
        self.N = N = bay.n_blocks

        self.b = b = model.binary_var_dict(product(irange(1, W), irange(1, H), irange(1, N), irange(1, N)), name='b')
        self.x = x = model.binary_var_dict(product(irange(1, W), irange(1, H), irange(1, W), irange(1, H), irange(1, N), irange(1, N)), name='x')
        self.y = y = model.binary_var_dict(product(irange(1, W), irange(1, H), irange(1, N), irange(1, N)), name='y')

        self.v = v = {(n, t): int(n < t) for n, t in product(irange(1, N), irange(1, N))}

        # objective
        model.minimize(model.sum(x[i, j, k, l, n, t] for i, j, k, l, n, t in product(irange(1, W), irange(1, H), irange(1, W), irange(1, H), irange(1, N), irange(1, N))))
        # (1)
        model.add_constraints(model.sum(b[i, j, n, t] for i, j in product(irange(1, W), irange(1, H))) + v[n, t] == 1 for n, t in product(irange(1, N), irange(1, N)))
        # (2)
        model.add_constraints(model.sum(b[i, j, n, t] for n in irange(1, N)) <= 1 for i, j, t in product(irange(1, W), irange(1, H), irange(1, N)))
        # (3)
        model.add_constraints(model.sum(b[i, j, n, t] for n in irange(1, N)) >= model.sum(b[i, j + 1, n, t] for n in irange(1, N)) for i, j, t in product(irange(1, W), irange(1, H - 1), irange(1, N)))
        # (6)
        model.add_constraints(b[i, j, n, t] == b[i, j, n, t - 1] + model.sum(x[k, l, i, j, n, t - 1] for k, l in product(irange(1, W), irange(1, H))) - model.sum(x[i, j, k, l, n, t - 1] for k, l in product(irange(1, W), irange(1, H))) - y[i, j, n, t - 1] for i, j, n, t in product(irange(1, W), irange(1, H), irange(1, N), irange(2, N)))
        # (7')
        model.add_constraints(v[n, t] == model.sum(y[i, j, n, tt] for i, j, tt in product(irange(1, W), irange(1, H), irange(1, t - 1))) for n, t in product(irange(1, N), irange(2, N)))
        # (8')
        model.add_constraints((H - 1) * (1 - model.sum(x[i, j, k, l, n, t] for n in irange(1, N))) >= model.sum(x[i, jj, k, ll, n, t] for n, jj, ll in product(irange(1, N), irange(j + 1, H), irange(l + 1, H))) for i, j, k, l, t in product(irange(1, W), irange(1, H - 1), irange(1, W), irange(1, H - 1), irange(1, N)))
        # (9)
        model.add_constraints((H - 1) * (1 - model.sum(b[i, j, t, t] for j in irange(1, H))) >= model.sum(x[ii, j, k, l, n, t] for ii in irange(1, W) if ii != i for j, k, l, n in product(irange(1, H), irange(1, W), irange(1, H), irange(1, N))) for i, t in product(irange(1, W), irange(1, N)))
        # (10)
        model.add_constraints(x[i, j, i, l, n, t] == 0 for i, j, l, n, t in product(irange(1, W), irange(1, H), irange(1, H), irange(1, N), irange(1, N)))
        # (A)
        model.add_constraints(model.sum(x[i, j + 1, k, l, n, t] for k, l, n in product(irange(1, W), irange(1, H), irange(1, N))) - model.sum(b[i, j + 1, n, t] for n in irange(1, N)) + 1 >= model.sum(x[i, j, k, l, n, t] for k, l, n in product(irange(1, W), irange(1, H), irange(1, N))) + model.sum(y[i, j, n, t] for n in irange(1, N)) for i, j, t in product(irange(1, W), irange(1, H - 1), irange(1, N)))
        # pre-processing
        model.add_constraints(b[i, j, n, 1] == int(bay.pri[i - 1][j - 1] == n) for i, j, n in product(irange(1, W), irange(1, H), irange(1, N)))

    def get_bays(self):
        bays = {}
        for t in irange(1, self.N):
            conf = [[None] * self.H for _ in range(self.W)]
            for i, j, n in product(irange(1, self.W), irange(1, self.H), irange(1, self.N)):
                if round(self.b[i, j, n, t].solution_value) == 1:
                    conf[i - 1][j - 1] = n
            bays[t] = Bay(self.W, self.H, conf)
        return bays

    def get_n_relos(self):
        return self.model.objective_value


def test():
    conf = [[2], [4, 6], [1, 3, 5]]
    bay = Bay(3, 3, conf)
    print('bay')
    print(bay)

    brp2c = BRP_II_C(bay)
    if brp2c.model.solve():
        print()
        print('n_relos = {}'.format(brp2c.get_n_relos()))
        bays = brp2c.get_bays()
        for t in irange(1, brp2c.N):
            print('t = {}'.format(t))
            print(bays[t])


if __name__ == '__main__':
    test()
