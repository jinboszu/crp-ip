`BRP_II_X.py` implements the BRP-II* model for the **restricted CRP with distinct priorities** in Expósito-Izquierdo et
al. (2015). The following bugs in BRP-II* are fixed:

- The index `t` in Constraints (10) should start with 2.
- Constraints (13) should be included.

Reference:

- Expósito-Izquierdo, C., Melián-Batista, B., & Moreno-Vega, J. M. (2015). An exact approach for the blocks relocation
  problem. *Expert Systems with Applications*, 42(17), 6408-6422.
