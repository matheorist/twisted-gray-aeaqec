# Twisted Gray Transfer and Singleton-Sharp Ebit--Dimension Frontiers

Reproducibility code for:

> Twisted Gray Transfer and Singleton-Sharp Ebit--Dimension Frontiers for Matrix-Product AEAQEC Codes from Split Semisimple F_q-Algebras

Author: Yong Zhang  
Affiliation: School of Mathematics and Statistics, Yancheng Teachers University, Yancheng 224002, P. R. China  
Email: zyyctc@126.com

## Purpose

These dependency-free Python scripts reproduce the finite-field matrix checks, two exact q=5, s=4 distance instances, the joint order/profile optimization over the complete q=5,n=4 budget range, the certified fixed-ebit frontier, the integer-envelope checks, and the inverse distance--ebit--dimension table reported in the manuscript. They are verification scripts rather than a general-purpose coding library.

## Requirements

Python 3.8 or later. No third-party packages are required. Run the scripts from this directory:

```text
python verify_AEAQEC.py
python verify_frontier_q5.py
python exact_distance_q5.py
python joint_optimize_q5.py
python verify_integer_envelope.py
python verify_inverse_region_q5.py
```

Each script prints `PASS` messages when its checks succeed.

## Version

The `v1.2.0` release contains the exact script versions used for the checks reported in the revised manuscript, including the complete joint-optimization curve and catalogue Pareto selection.

## License

MIT License; see [LICENSE](LICENSE).
