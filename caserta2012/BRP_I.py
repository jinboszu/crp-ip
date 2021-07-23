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


class BRP_I:
    def __init__(self, bay):
        bay.validate_full_distinct()
        self.model = model = Model()

        self.W = W = bay.n_stacks
        self.H = H = bay.n_tiers
        self.N = N = bay.n_blocks
        self.T = T = N + bay.brp_min_max()

        self.b = b = model.binary_var_dict(product(irange(1, W), irange(1, H), irange(1, N), irange(1, T)), name='b')
        self.v = v = model.binary_var_dict(product(irange(1, N), irange(1, T)), name='v')
        self.x = x = model.binary_var_dict(product(irange(1, W), irange(1, H), irange(1, W), irange(1, H), irange(1, N), irange(1, T)), name='x')
        self.y = y = model.binary_var_dict(product(irange(1, W), irange(1, H), irange(1, N), irange(1, T)), name='y')

        # objective
        model.maximize(model.sum(v[N, t] for t in irange(1, T)))
        # (1)
        model.add_constraints(model.sum(b[i, j, n, t] for i, j in product(irange(1, W), irange(1, H))) + v[n, t] == 1 for n, t in product(irange(1, N), irange(1, T)))
        # (2)
        model.add_constraints(model.sum(b[i, j, n, t] for n in irange(1, N)) <= 1 for i, j, t in product(irange(1, W), irange(1, H), irange(1, T)))
        # (3)
        model.add_constraints(model.sum(b[i, j, n, t] for n in irange(1, N)) >= model.sum(b[i, j + 1, n, t] for n in irange(1, N)) for i, j, t in product(irange(1, W), irange(1, H - 1), irange(1, T)))
        # (4)
        model.add_constraints(model.sum(x[i, j, k, l, n, t] for i, j, k, l, n in product(irange(1, W), irange(1, H), irange(1, W), irange(1, H), irange(1, N))) + model.sum(y[i, j, n, t] for i, j, n in product(irange(1, W), irange(1, H), irange(1, N))) <= 1 for t in irange(1, T))
        # (5)
        model.add_constraints(model.sum(v[n, t] for t in irange(1, T)) >= model.sum(v[n + 1, t] for t in irange(1, T)) + 1 for n in irange(1, N - 1))
        # (6)
        model.add_constraints(b[i, j, n, t] == b[i, j, n, t - 1] + model.sum(x[k, l, i, j, n, t - 1] for k, l in product(irange(1, W), irange(1, H))) - model.sum(x[i, j, k, l, n, t - 1] for k, l in product(irange(1, W), irange(1, H))) - y[i, j, n, t - 1] for i, j, n, t in product(irange(1, W), irange(1, H), irange(1, N), irange(2, T)))
        # (7)
        model.add_constraints(v[n, t] == model.sum(y[i, j, n, tt] for i, j, tt in product(irange(1, W), irange(1, H), irange(1, t - 1))) for n, t in product(irange(1, N), irange(1, T)))
        # pre-processing
        model.add_constraints(b[i, j, n, 1] == int(bay.pri[i - 1][j - 1] == n) for i, j, n in product(irange(1, W), irange(1, H), irange(1, N)))

    def get_bays(self):
        bays = {}
        for t in irange(1, self.T):
            conf = [[None] * self.H for _ in range(self.W)]
            for i, j, n in product(irange(1, self.W), irange(1, self.H), irange(1, self.N)):
                if round(self.b[i, j, n, t].solution_value) == 1:
                    conf[i - 1][j - 1] = n
            bays[t] = Bay(self.W, self.H, conf)
        return bays

    def get_n_relos(self):
        return sum(round(self.x[i, j, k, l, n, t].solution_value) for i, j, k, l, n, t in product(irange(1, self.W), irange(1, self.H), irange(1, self.W), irange(1, self.H), irange(1, self.N), irange(1, self.T)))


def test():
    conf = [[1, 3, 4], [5], [2]]
    bay = Bay(3, 3, conf)
    print('bay')
    print(bay)

    brp_i = BRP_I(bay)
    if brp_i.model.solve():
        print()
        print('n_relos = {}'.format(brp_i.get_n_relos()))
        bays = brp_i.get_bays()
        for t in irange(1, brp_i.T):
            print('t = {}'.format(t))
            print(bays[t])


if __name__ == '__main__':
    test()
