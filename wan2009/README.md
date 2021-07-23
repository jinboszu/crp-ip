`MRIP.py` implements the MRIP model for the **restricted CRP with distinct priorities** in Wan et al. (2009). The
following bug in MRIP is fixed:

- The index `i` starts with 2 in Constraints (26) in Wan et al. (2009). In such case, `x_{11cp}` are free, which may be
  different from `X_{11cp}` in the solution obtained. In my implementation, Constraints (26) is modified to ensure
  `x_{11cp} = X_{11cp}`.

Reference:

- Wan, Y.-w., Liu, J., & Tsai, P.-C. (2009). The assignment of storage locations to containers for a container stack.
  *Naval Research Logistics (NRL)*, 56(8), 699-713.
