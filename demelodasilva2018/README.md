`BRP_m1.py` and `BRP_m2.py` implement the BRP-m1 and BRP-m2 models for the **unrestricted CRP with duplicate
priorities**, respectively, in de Melo da Silva et al. (2018). The following bugs in BRP-m2 are fixed:

- Retrievals are not carried out in time 0. BRP-m2 may fail when there are retrievable blocks at the beginning and `T`
  is tightly set.
- The right-hand-side in Constraints (49) (when `G < N`) should be `1-\sum_{l=1}^{g-1}(x^t_{lsi}+k^t_{lsi})`, otherwise
  a stack `[1, 2, 3]` (bottom up) can be immediately cleared.
- Relocations and retrievals are not forced to be carried out as early as possible. Two groups of constraints are added
  for this purpose.

Note that the **restricted** versions of BRP-m1 and BRP-m2, referred to as R-BRP-m1 and R-BRP-m2, respectively, are only
applicable to **distinct priorities**. The following bug in R-BRP-m2 is fixed:

- The first term of the right-hand-side in Constriants (52) has to represent whether block `g` has been retrieved
  before (not after) time `t`, and the second term has to represent whether a block `i <= g` exists in stack `s` before
  (not after) time `t`.

Reference:

- de Melo da Silva, M., Toulouse, S., & Wolfler Calvo, R. (2018). A new effective unified model for solving the
  pre-marshalling and block relocation problems. *European Journal of Operational Research*, 271(1), 40-56.
