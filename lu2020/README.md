`BRP_m3.py` implements the BRP-m3 model for the **unrestricted CRP with duplicate priorities** in Lu et al. (2020). The
following bug in BRP-m3 is fixed:

- Retrievals are not carried out in time 0. BRP-m2 may fail when there are retrievable blocks at the beginning and `T`
  is tightly set.

Note that Constraints (e1)-(e7) have been modified to apply to duplicate priorities.

Reference:

- Lu, C., Zeng, B., & Liu, S. (2020). A study on the block relocation problem: Lower bound derivations and strong
  formulations. *IEEE Transactions on Automation Science and Engineering*, 17(4), 1829-1853.
