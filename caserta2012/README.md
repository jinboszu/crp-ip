`BRP_I.py` and `BRP_II.py` implement the BRP-I and BRP-II models for the **unrestricted and restricted CRP with distinct
priorities**, respectively, in Caserta et al. (2012). The following bugs in BRP-II are fixed:

- Constraints (8) in BRP-II do not allow relocating several blocks to the same destination stack.
- Constraints (9) in BRP-II allow relocating blocks below the target block.

Reference:

- Caserta, M., Schwarze, S., & Voß, S. (2012). A mathematical formulation and complexity considerations for the blocks
  relocation problem. *European Journal of Operational Research*, 219(1), 96-104.
