# Twisted Gray Transfer and Singleton-Sharp Ebit--Dimension Frontiers

Reproducibility code for:

> Twisted Gray Transfer and Singleton-Sharp Ebit--Dimension Frontiers for Matrix-Product AEAQEC Codes from Split Semisimple F_q-Algebras

Author: Yong Zhang  
Affiliation: School of Mathematics and Statistics, Yancheng Teachers University, Yancheng 224002, P. R. China  
Email: zyyctc@126.com

## Purpose

These dependency-free Python scripts reproduce the finite-field matrix checks, two exact q=5, s=4 distance instances, the full three-point antichain-boundary lift, and a five-point exact-distance segment of the K=10 plateau, the joint order/profile optimization over the complete q=5,n=4 budget range, the target-aware distance-grid search, its matrix-independent separate-cost prefilter, the certified fixed-ebit frontier, the integer-envelope checks, and the inverse distance--ebit--dimension table reported in the manuscript. They are verification scripts rather than a general-purpose coding library.

## Requirements

Python 3.8 or later. No third-party packages are required. Run the scripts from this directory:

```text
python verify_AEAQEC.py
python verify_frontier_q5.py
python exact_distance_q5.py
python verify_witness_coverage_q5.py
python joint_optimize_q5.py
python joint_target_optimize_q5.py
python verify_integer_envelope.py
python verify_inverse_region_q5.py
python verify_inverse_budget_q5.py
```

Each script prints `PASS` messages when its checks succeed. The witness-coverage script verifies exact (d_Z,d_X)=(3,3) at c=0,1,2,3,4 with K=10 and exact (4,2) and (2,4) boundary witnesses at c=0, so the complete target antichain {(4,2),(3,3),(2,4)} is covered at K=10. The target-aware script also verifies the threshold 12 prefilter, the downward-closed target region, and its maximal antichain; it prunes 34 of the 64 pairs in the 1..8 target grid before any matrix search. The inverse-budget script computes the exact minimum ebit count for every 1..4 target pair at K=9 and checks the separate-cost lower bound and the gate ceiling/deficit for the unavailable (4,4) target and the three-point antichain boundary; it also verifies the K=9 boundary transition at cmax=0,1,12.

## Version

The `v2.0.0` release contains the exact script versions used for the checks reported in the revised manuscript, including the complete antichain-boundary exact witnesses and the five-point exact-distance witness-coverage segment, complete joint-optimization curve, target-aware search with the matrix-independent prefilter, exact inverse-budget table, gate-deficit certificate, target-antichain boundary, budget activation events, and catalogue Pareto selection.

## License

MIT License; see [LICENSE](LICENSE).











