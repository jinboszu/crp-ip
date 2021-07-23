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


from copy import copy
from itertools import product, groupby

from docplex.mp.model import Model

from bay import Bay
from common import irange


class BRP_m3:
    def __init__(self, bay, bug_fixed=True, restricted=False):
        self.model = model = Model()

        self.S = S = bay.n_stacks
        self.H = H = bay.n_tiers
        self.B = B = bay.n_blocks
        self.T = T = bay.brp_min_max()
        self.stack = stack = dict(zip(irange(1, B), (s + 1 for s in range(bay.n_stacks) for _ in range(bay.h[s]))))
        self.tier = tier = dict(zip(irange(1, B), (t + 1 for s in range(bay.n_stacks) for t in range(bay.h[s]))))
        self.p = p = dict(zip(irange(1, B), (bay.pri[s][t] for s in range(bay.n_stacks) for t in range(bay.h[s]))))
        p[B + 1] = bay.p_max + 1

        self.x = x = model.continuous_var_dict(((t, i, j) for t in irange(0, T) for i, j in product(irange(1, B), irange(1, B + 1)) if j != i), lb=0, ub=1, name='x')
        self.yout = yout = model.binary_var_dict(((t, i, j) for t in irange(1, T) for i, j in product(irange(1, B), irange(1, B + 1)) if j != i), name='yout')
        self.yin = yin = model.binary_var_dict(((t, i, j) for t in irange(1, T) for i, j in product(irange(1, B), irange(1, B + 1)) if j != i), name='yin')
        self.z = z = model.binary_var_dict(((t, i, j) for t in irange(0, T) for i in irange(1, B) for j in irange(1, B + 1) if j != i and p[j] >= p[i]), name='z')
        self.u = u = model.continuous_var_dict(product(irange(1, T), irange(1, B)), lb=1, ub=H, name='u')

        G = {v: list(l) for v, l in groupby(sorted(irange(1, B), key=lambda i: p[i]), key=lambda i: p[i])}
        C = {(i, j): int(stack[i] == stack[j] and tier[i] == tier[j] + 1 if j <= B else int(tier[i] == 1)) for i, j in product(irange(1, B), irange(1, B + 1)) if j != i}
        L = bay.compute_lb_kh()

        # objective
        model.minimize(model.sum(yin[t, i, j] for i, j in product(irange(1, B), irange(1, B + 1)) if j != i for t in irange(1, T)))
        # (x1)
        if bug_fixed:
            model.add_constraints(x[0, i, j] + (z[0, i, j] if p[j] >= p[i] else 0) == C[i, j] for i, j in product(irange(1, B), irange(1, B + 1)) if j != i)
        else:
            model.add_constraints(x[0, i, j] == C[i, j] for i, j in product(irange(1, B), irange(1, B + 1)) if j != i)
            model.add_constraints(z[0, i, j] == 0 for i, j in product(irange(1, B), irange(1, B + 1)) if j != i)
        # (x2*)
        model.add_constraints(x[t, i, j] == x[t - 1, i, j] - yout[t, i, j] + yin[t, i, j] for i, j in product(irange(1, B), irange(1, B + 1)) if p[j] < p[i] for t in irange(1, T))
        # (x3*)
        model.add_constraints(x[t, i, j] == x[t - 1, i, j] - yout[t, i, j] + yin[t, i, j] - z[t, i, j] for i, j in product(irange(1, B), irange(1, B + 1)) if j != i and p[j] >= p[i] for t in irange(1, T))
        # (x4)
        model.add_constraints(x[T, i, j] == 0 for i, j in product(irange(1, B), irange(1, B + 1)) if j != i)
        # (yout1)
        model.add_constraints(model.sum(yout[t, i, j] for i, j in product(irange(1, B), irange(1, B + 1)) if j != i) == 1 for t in irange(1, L))
        # (yout2)
        model.add_constraints(model.sum(yout[t, i, j] for i, j in product(irange(1, B), irange(1, B + 1)) if j != i) <= model.sum(yout[t - 1, i, j] for i, j in product(irange(1, B), irange(1, B + 1)) if j != i) for t in irange(max(2, L + 1), T))
        # (yout4)
        model.add_constraints(model.sum(yout[t, i, j] for j in irange(1, B + 1) if j != i) <= model.sum(x[t - 1, i, j] for j in irange(1, B + 1) if j != i) - model.sum(x[t - 1, j, i] for j in irange(1, B) if j != i) for i, t in product(irange(1, B), irange(1, T)))
        # (yin1)
        model.add_constraints(model.sum(yin[t, i, j] for j in irange(1, B + 1) if j != i) == model.sum(yout[t, i, j] for j in irange(1, B + 1) if j != i) for i, t in product(irange(1, B), irange(1, T)))
        # (yin2)
        model.add_constraints(model.sum(yin[t, j, i] for j in irange(1, B) if j != i) <= model.sum(x[t - 1, i, j] for j in irange(1, B + 1) if j != i) - model.sum(yout[t, j, i] for j in irange(1, B) if j != i) for i, t in product(irange(1, B), irange(1, T)))
        # (yin3)
        model.add_constraints(model.sum(yin[t, j, B + 1] for j in irange(1, B)) <= 1 - model.sum(yout[t, j, B + 1] for j in irange(1, B)) for t in irange(1, T))
        # (yin4)
        model.add_constraints(model.sum(yin[t, j, i] for j in irange(1, B) if j != i) <= model.sum(x[t - 1, i, j] for j in irange(1, B + 1) if j != i) - model.sum(x[t - 1, j, i] for j in irange(1, B) if j != i) for i, t in product(irange(1, B), irange(1, T)))
        # (yin5)
        model.add_constraints(model.sum(yin[t, j, B + 1] for j in irange(1, B)) <= S - model.sum(x[t - 1, j, B + 1] for j in irange(1, B)) for t in irange(1, T))
        # (z3*)
        if bug_fixed:
            model.add_constraints(len(G.get(p[i] - 1, [])) * model.sum(z[tt, i, j] for j in irange(1, B + 1) if j != i and p[j] >= p[i] for tt in irange(0, t)) <= model.sum(z[tt, k, j] for k in G.get(p[i] - 1, []) for j in irange(1, B + 1) if j != k and p[j] >= p[k] for tt in irange(0, t)) for i in irange(1, B) if p[i] >= 2 for t in irange(0, T))
        else:
            model.add_constraints(len(G.get(p[i] - 1, [])) * model.sum(z[tt, i, j] for j in irange(1, B + 1) if j != i and p[j] >= p[i] for tt in irange(1, t)) <= model.sum(z[tt, k, j] for k in G.get(p[i] - 1, []) for j in irange(1, B + 1) if j != k and p[j] >= p[k] for tt in irange(1, t)) for i in irange(1, B) if p[i] >= 2 for t in irange(1, T))
        # (u1)
        model.add_constraints(u[t, i] >= u[t, j] + 1 - H * (1 - x[t - 1, i, j] + yout[t, i, j] - yin[t, i, j]) for i, j in product(irange(1, B), irange(1, B)) if j != i for t in irange(1, T))
        # (e1*)
        model.add_constraints(z[t, k, i] == (yout[t, j, k] if t > 0 else 0) + (z[t, j, k] if p[j] == 1 else 0) for k in G.get(1, []) for i in irange(1, B + 1) if i != k and C[k, i] == 1 for j in irange(1, B) if j != k and C[j, k] == 1 for t in irange(0, T))
        # (e2*)
        model.add_constraints(model.sum((yout[t, i, k] if t > 0 else 0) + (z[t, i, k] if p[i] == 1 else 0) for t in irange(0, T)) == 1 for k in G.get(1, []) for i in irange(1, B) if i != k and C[i, k] == 1)
        # (e3*)
        model.add_constraints(yin[t, i, k] == 0 for k in G.get(1, []) for i in irange(1, B) if i != k for t in irange(1, T))
        # (e4*)
        model.add_constraints(x[t, i, k] == 0 for k in G.get(1, []) for i in irange(1, B) if i != k and C[i, k] == 0 for t in irange(0, T))
        model.add_constraints(yout[t, i, k] == 0 for k in G.get(1, []) for i in irange(1, B) if i != k and C[i, k] == 0 for t in irange(1, T))
        # (e5*)
        model.add_constraints(yout[t, k, i] == 0 for k in G.get(1, []) for i in irange(1, B + 1) if i != k for t in irange(1, T))
        model.add_constraints(yin[t, k, i] == 0 for k in G.get(1, []) for i in irange(1, B + 1) if i != k for t in irange(1, T))
        # (e6*)
        model.add_constraints(x[t, k, i] == 0 for k in G.get(1, []) for i in irange(1, B + 1) if i != k and C[k, i] == 0 for t in irange(0, T))
        model.add_constraints(z[t, k, i] == 0 for k in G.get(1, []) for i in irange(1, B + 1) if i != k and C[k, i] == 0 for t in irange(0, T))
        # (e7*)
        model.add_constraints(u[t, k] == tier[k] for k, t in product(G.get(1, []), irange(1, T)))

        if restricted:
            # (yout6)
            model.add(model.sum(yout[t, i, j] for j in irange(1, B + 1) if j != i) >= model.sum(yout[t - 1, j, i] for j in irange(1, B) if j != i) - model.sum(z[t - 1, i, j] for j in irange(1, B + 1) if j != i and p[j] >= p[i]) for i, t in product(irange(1, B), irange(2, T)))
            # (yout7)
            model.add_constraints(yout[t, t, B + 1] == 0 for i, t in product(irange(1, B), irange(1, T)))

    def get_bays(self):
        stack = copy(self.stack)
        tier = copy(self.tier)
        h = {s: 0 for s in irange(1, self.S)}
        for i in irange(1, self.B):
            h[stack[i]] += 1

        before = {}
        after = {}
        for t in irange(0, self.T):
            if t != 0:
                for i, j in product(irange(1, self.B), irange(1, self.B + 1)):
                    if j != i and round(self.yin[t, i, j].solution_value) == 1:
                        if j <= self.B:
                            stack[i] = stack[j]
                            tier[i] = tier[j] + 1
                        else:
                            stack[i] = next(s for s in irange(1, self.S) if h[s] == 0)
                            tier[i] = 1
            conf_before = [[None] * self.H for _ in range(self.S)]
            for i in irange(1, self.B):
                if sum(round(self.x[t, i, j].solution_value) + (round(self.z[t, i, j].solution_value) if self.p[j] >= self.p[i] else 0) for j in irange(1, self.B + 1) if j != i) == 1:
                    conf_before[stack[i] - 1][tier[i] - 1] = self.p[i]

            conf_after = [[None] * self.H for _ in range(self.S)]
            for i in irange(1, self.B):
                if sum(round(self.x[t, i, j].solution_value) for j in irange(1, self.B + 1) if j != i) == 1:
                    conf_after[stack[i] - 1][tier[i] - 1] = self.p[i]

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

    brp_m3 = BRP_m3(bay)
    if brp_m3.model.solve():
        print()
        print('n_relos = {}'.format(brp_m3.get_n_relos()))
        before, after = brp_m3.get_bays()
        for t in irange(0, brp_m3.T):
            print('t = {} (before retrievals)'.format(t))
            print(before[t])
            print('t = {} (after retrievals)'.format(t))
            print(after[t])


if __name__ == '__main__':
    test()
