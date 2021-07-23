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


from copy import copy, deepcopy
from itertools import chain, islice

from common import irange


class Bay:
    def __init__(self, n_stacks, n_tiers, conf):
        assert len(conf) == n_stacks
        for s in range(n_stacks):
            assert len(conf[s]) <= n_tiers
            for t in range(1, len(conf[s])):
                assert conf[s][t] is None or conf[s][t - 1] is not None

        self.n_stacks = n_stacks
        self.n_tiers = n_tiers
        self.h = h = [len(conf[i]) - conf[i].count(None) for i in range(n_stacks)]
        self.n_blocks = sum(self.h)
        self.pri = pri = [conf[s] + [None] * (n_tiers - len(conf[s])) for s in range(n_stacks)]
        self.values = values = list(chain(*(islice(pri[s], h[s]) for s in range(n_stacks))))
        self.p_max = p_max = max(values, default=0)
        self.qlt = qlt = [[None] * n_tiers + [p_max + 1] for _ in range(n_stacks)]

        for s in range(n_stacks):
            for t in range(h[s]):
                qlt[s][t] = min(qlt[s][t - 1], pri[s][t])

    def is_distinct(self):
        values = list(chain(*(islice(self.pri, self.h[s]) for s in range(self.n_stacks))))
        return len(values) == 0 or len(values) == len(set(values)) and min(values) + len(values) - 1 == max(values)

    def __str__(self):
        n_stacks = self.n_stacks
        n_tiers = self.n_tiers
        h = self.h
        pri = self.pri
        p_max = self.p_max

        tier_width = len(str(n_tiers))
        cell_width = len(str(max(n_stacks, p_max)))

        builder = []

        h_max = max(h)
        if h_max <= n_tiers - 1:
            builder.append(str(n_tiers).rjust(tier_width))
            builder.append(' | ')
            builder.append('[{}]'.format(' ' * cell_width) * n_stacks)
            builder.append('\n')

        if h_max in [n_tiers - 2, n_tiers - 3]:
            builder.append(str(n_tiers - 1).rjust(tier_width))
            builder.append(' | ')
            builder.append('[{}]'.format(' ' * cell_width) * n_stacks)
            builder.append('\n')

        if h_max <= n_tiers - 4:
            builder.append('.' * tier_width)
            builder.append(' | ')
            builder.append(' {} '.format('.' * cell_width) * n_stacks)
            builder.append('\n')

        if h_max <= n_tiers - 3:
            builder.append(str(h_max + 1).rjust(tier_width))
            builder.append(' | ')
            builder.append('[{}]'.format(' ' * cell_width) * self.n_stacks)
            builder.append('\n')

        for t in reversed(range(h_max)):
            builder.append(str(t + 1).rjust(tier_width))
            builder.append(' | ')
            builder.append(''.join('[{}]'.format(str(pri[s][t]).rjust(cell_width)) if t < h[s] else '[{}]'.format(' ' * cell_width) for s in range(n_stacks)))
            builder.append('\n')

        builder.append('-' * tier_width)
        builder.append('-' * 3)
        builder.append('-' * ((cell_width + 2) * n_stacks))
        builder.append('\n')

        builder.append(' ' * tier_width)
        builder.append(' | ')
        builder.append(''.join(' {} '.format(str(s + 1).rjust(cell_width)) for s in range(n_stacks)))

        return ''.join(builder)

    def brp_min_max(self):
        if self.n_blocks == 0:
            return 0

        n_stacks = self.n_stacks
        n_tiers = self.n_tiers
        n_blocks = self.n_blocks
        h = copy(self.h)
        pri = deepcopy(self.pri)
        qlt = deepcopy(self.qlt)

        n_relos = 0
        while n_blocks > 0:
            min_value = min(chain(*(islice(pri[s], h[s]) for s in range(n_stacks))))
            (s_target, t_target) = next((s, t) for s in range(n_stacks) for t in reversed(range(h[s])) if pri[s][t] == min_value and t + 1 >= n_blocks - (n_stacks - 1) * n_tiers)
            while t_target < h[s_target] - 1:
                p = pri[s_target][h[s_target] - 1]

                s_dest = min((s for s in range(n_stacks) if s != s_target and h[s] < n_tiers and p <= qlt[s][h[s] - 1]), key=lambda s: qlt[s][h[s] - 1], default=max((s for s in range(n_stacks) if s != s_target and h[s] < n_tiers), key=lambda s: qlt[s][h[s] - 1]))

                pri[s_target][h[s_target] - 1] = None
                qlt[s_target][h[s_target] - 1] = None
                h[s_target] -= 1

                h[s_dest] += 1
                pri[s_dest][h[s_dest] - 1] = p
                qlt[s_dest][h[s_dest] - 1] = min(qlt[s_dest][h[s_dest] - 2], p)

                n_relos += 1

            pri[s_target][h[s_target] - 1] = None
            qlt[s_target][h[s_target] - 1] = None
            h[s_target] -= 1
            n_blocks -= 1

        return n_relos

    def validate_full_distinct(self):
        values = sorted(self.values)
        for i in range(self.n_blocks):
            assert values[i] == i + 1

    def validate_distinct(self):
        assert len(set(self.values)) == self.n_blocks

    def compute_lb_kh(self):
        n_bad = 0
        for s in range(self.n_stacks):
            for t in range(self.h[s]):
                n_bad += int(self.pri[s][t] > self.qlt[s][t])
        return n_bad

    def compute_lb_zhu(self):
        self.validate_distinct()
        n_stacks = self.n_stacks
        n_tiers = self.n_tiers
        n_blocks = self.n_blocks
        pri = self.pri
        qlt = self.qlt

        h = copy(self.h)
        loc = {pri[s][t]: (s, t) for s in range(n_stacks) for t in range(h[s])}
        q_max = max(qlt[s][h[s] - 1] for s in range(n_stacks) if h[s] < n_tiers)

        lb = {}
        lb_plus = {}
        for n in irange(1, n_blocks):
            (s, t) = loc[n]
            if t < h[s]:
                lb[n] = h[s] - t - 1
                lb_plus[n] = sum(1 for tt in range(t + 1, h[s]) if pri[s][tt] > q_max)
                h[s] = t
                q_max = max(q_max, qlt[s][h[s] - 1])
            else:
                lb[n] = 0
                lb_plus[n] = 0

        return lb, lb_plus
