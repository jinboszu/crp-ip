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


class BRP_II_A:
    def __init__(self, bay):
        bay.validate_full_distinct()
        self.model = model = Model()

        self.W = W = bay.n_stacks
        self.H = H = bay.n_tiers
        self.N = N = bay.n_blocks

        self.b = b = model.binary_var_dict(((i, j, n, t) for i, j, t in product(irange(1, W), irange(1, H), irange(1, N)) for n in irange(t, N)), name='b')
        self.x = x = model.binary_var_dict(((i, j, k, l, n, t) for i, j, k, l, t in product(irange(1, W), irange(2, H), irange(1, W), irange(1, H), irange(1, N - 1)) for n in irange(t + 1, N)), name='x')
        self.y = y = model.binary_var_dict(((i, j, t, t) for i, j, t in product(irange(1, W), irange(1, H), irange(1, N - 1))), name='y')

        lb, lb_plus = bay.compute_lb_zhu()
        lb_sum = model.sum(lb.values())
        lb_plus_right = {}
        for t in reversed(irange(1, N)):
            lb_plus_right[t] = (lb_plus_right[t + 1] if t < N else 0) + lb_plus[t]
        ub = bay.brp_min_max()
        UB = {t: min(ub - 1 - lb_sum + lb[t] - lb_plus_right[t], H - 1) for t in irange(1, N)}
        Q = {bay.pri[s][t]: bay.qlt[s][t] for s in range(bay.n_stacks) for t in range(bay.h[s])}
        ix = {bay.pri[s][t]: s + 1 for s in range(bay.n_stacks) for t in range(bay.h[s])}
        jx = {bay.pri[s][t]: t + 1 for s in range(bay.n_stacks) for t in range(bay.h[s])}

        # objective
        model.minimize(model.sum(x[i, j, k, l, n, t] for i, j, k, l, t in product(irange(1, W), irange(2, H), irange(1, W), irange(1, H), irange(1, N - 1)) for n in irange(t + 1, N)))
        # (2)
        model.add_constraints(model.sum(b[i, j, n, t] for n in irange(t, N)) <= 1 for i, j, t in product(irange(1, W), irange(1, H), irange(1, N - 1)))
        # (3)
        model.add_constraints(model.sum(b[i, j, n, t] for n in irange(t, N)) >= model.sum(b[i, j + 1, n, t] for n in irange(t, N)) for i, j, t in product(irange(1, W), irange(1, H - 1), irange(1, N)))
        # (6a)
        model.add_constraints(b[i, j, n, t + 1] == b[i, j, n, t] + model.sum(x[k, l, i, j, n, t] for k, l in product(irange(1, W), irange(2, H))) - (0 if j == 1 else model.sum(x[i, j, k, l, n, t] for k, l in product(irange(1, W), irange(1, H)))) for i, j, t in product(irange(1, W), irange(1, H), irange(1, N - 1)) for n in irange(t + 1, N))
        # (6b)
        model.add_constraints(b[i, j, t, t] - y[i, j, t, t] == 0 for i, j, t in product(irange(1, W), irange(1, H), irange(1, N - 1)))
        # (7'')
        model.add_constraints(model.sum(y[i, j, t, t] for i, j in product(irange(1, W), irange(1, H))) == 1 for t in irange(1, N - 1))
        # (8')
        model.add_constraints((H - 1) * (1 - model.sum(x[i, j, k, l, n, t] for n in irange(t + 1, N))) >= model.sum(x[i, jj, k, ll, n, t] for n, jj, ll in product(irange(t + 1, N), irange(j + 1, H), irange(l + 1, H))) for i, j, k, l, t in product(irange(1, W), irange(2, H - 1), irange(1, W), irange(1, H - 1), irange(1, N - 1)))
        # (A')
        model.add_constraints(model.sum(y[i, jj, t, t] for jj in irange(1, j - 1)) >= model.sum(x[i, j, k, l, n, t] for k, l, n in product(irange(1, W), irange(1, H), irange(t + 1, N))) for i, j, t in product(irange(1, W), irange(2, H), irange(1, N - 1)))
        # (B)
        model.add_constraints(model.sum(x[i, j, k, l, n, t] for i, j, k, l, n in product(irange(1, W), irange(2, H), irange(1, W), irange(1, H), irange(t + 1, N))) <= UB[t] for t in irange(1, N - 1))
        # pre-processing
        # Block n is at position (i_n, j_n) and nowhere else until period π_n.
        # No other block n' may be located at position (i_n, j_n) until period π_n.
        model.add_constraints(b[i, j, n, t] == int(bay.pri[i - 1][j - 1] == n) for i, j, n in product(irange(1, W), irange(1, H), irange(1, N)) for t in irange(1, Q[n]))
        # Position (i_n, j_n) is occupied by block n until period π_n.
        # Hence, no other block n' < n may be retrieved from position (i_n, j_n) prior to period π_n.
        # No other block n' may be relocated from or to position (i_n, j_n) until period π_n.
        model.add_constraints(y[ix[n], jx[n], t, t] == 0 for n in irange(1, N) for t in irange(1, Q[n] - 1))
        model.add_constraints(x[ix[n], jx[n], k, l, nn, t] == 0 for n in irange(1, N) if jx[n] >= 2 for nn in irange(2, N) if nn != n for k, l, t in product(irange(1, W), irange(1, H), irange(1, min(nn - 1, Q[n]))))
        model.add_constraints(x[i, j, ix[n], jx[n], nn, t] == 0 for n in irange(1, N) for nn in irange(2, N) if nn != n for i, j, t in product(irange(1, W), irange(2, H), irange(1, min(nn - 1, Q[n]))))
        # If π_n < n, block n is relocated for the first time in period π_n.
        # It is not relocated prior to period π_n or from a position other than (i_n, j_n) in period π_n.
        model.add_constraints(x[i, j, k, l, n, t] == 0 for n in irange(2, N) if Q[n] >= 2 for i, j, k, l, t in product(irange(1, W), irange(2, H), irange(1, W), irange(1, H), irange(1, min(n - 1, Q[n] - 1))))
        model.add_constraints(x[i, j, k, l, n, Q[n]] == 0 for n in irange(2, N) if Q[n] < n for i, j in product(irange(1, W), irange(2, H)) if bay.pri[i - 1][j - 1] != n for k, l in product(irange(1, W), irange(1, H)))
        # If π_n = n, block n is never relocated and retrieved from its initial position (i_n, j_n) in period t = n.
        model.add_constraints(x[i, j, k, l, n, t] == 0 for n in irange(2, N) if Q[n] == n for i, j, k, l, t in product(irange(1, W), irange(2, H), irange(1, W), irange(1, H), irange(1, n - 1)))
        model.add_constraints(y[i, j, t, t] == 0 for t in irange(1, N - 1) if Q[t] == t for i, j in product(irange(1, W), irange(1, H)) if bay.pri[i - 1][j - 1] != t)
        # If π_n = n, only blocks in stack i_n and above j_n may be relocated in period π_n (Assumption A1).
        model.add_constraints(x[i, j, k, l, nn, t] == 0 for t in irange(1, N - 1) if Q[t] == t for i in irange(1, W) if i != ix[t] for j, k, l, nn in product(irange(2, H), irange(1, W), irange(1, H), irange(t + 1, N)))
        model.add_constraints(x[ix[t], j, k, l, nn, t] == 0 for t in irange(1, N - 1) if Q[t] == t and jx[t] >= 2 for j, k, l, nn in product(irange(2, jx[t]), irange(1, W), irange(1, H), irange(t + 1, N)))
        # If π_n = n, position (i_n, j_n) and all positions above have to be empty in period n + 1.
        # No containers can be retrieved or relocated from these positions in period π_n + 1.
        model.add_constraints(b[ix[n], j, nn, n + 1] == 0 for n in irange(1, N - 1) if Q[n] == n for j, nn in product(irange(max(2, jx[n]), H), irange(n + 1, N)))
        model.add_constraints(y[ix[n], j, n + 1, n + 1] == 0 for n in irange(1, N - 2) if Q[n] == n for j in irange(max(2, jx[n]), H))
        model.add_constraints(x[ix[n], j, k, l, nn, n + 1] == 0 for n in irange(1, N - 2) if Q[n] == n for j, k, l, nn in product(irange(max(2, jx[n]), H), irange(1, W), irange(1, H), irange(n + 2, N)))
        # At period t, tiers at height h > N_t may not be occupied and no containers can be retrieved or relocated from these positions:
        model.add_constraints(b[i, j, n, t] == 0 for t in irange(max(1, N + 2 - H), N) for i, j, n in product(irange(1, W), irange(N + 2 - t, H), irange(t, N)))
        model.add_constraints(y[i, j, t, t] == 0 for t in irange(max(1, N + 2 - H), N - 1) for i, j in product(irange(1, W), irange(N + 2 - t, H)))
        model.add_constraints(x[i, j, k, l, n, t] == 0 for t in irange(max(1, N + 2 - H), N - 1) for i, j, k, l, n in product(irange(1, W), irange(N + 2 - t, H), irange(1, W), irange(1, H), irange(t + 1, N)))
        # At period t, relocation blocks can only be put into tiers h <= N_{t + 1} = N_t − 1:
        model.add_constraints(x[i, j, k, l, n, t] == 0 for t in irange(max(1, N + 1 - H), N - 1) for i, j, k, l, n in product(irange(1, W), irange(2, H), irange(1, W), irange(N + 1 - t, H), irange(t + 1, N)))
        # Relocations x_{ijklnt} with i = k may not exist.
        model.add_constraints(x[i, j, i, l, n, t] == 0 for i, j, l, n in product(irange(1, W), irange(2, H), irange(1, W), irange(2, N)) for t in irange(1, n - 1))

    def get_bays(self):
        bays = {}
        for t in irange(1, self.N):
            mat = [[None] * self.H for _ in range(self.W)]
            for i, j, n in product(irange(1, self.W), irange(1, self.H), irange(t, self.N)):
                if round(self.b[i, j, n, t].solution_value) == 1:
                    mat[i - 1][j - 1] = n
            bays[t] = Bay(self.W, self.H, mat)
        return bays

    def get_n_relos(self):
        return self.model.objective_value


def test():
    conf = [[2], [4, 6], [1, 3, 5]]
    bay = Bay(3, 3, conf)
    print('bay')
    print(bay)

    brp2c = BRP_II_A(bay)
    if brp2c.model.solve():
        print()
        print('n_relos = {}'.format(brp2c.get_n_relos()))
        bays = brp2c.get_bays()
        for t in irange(1, brp2c.N):
            print('t = {}'.format(t))
            print(bays[t])


if __name__ == '__main__':
    test()
