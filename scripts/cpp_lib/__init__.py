"""
cpp_lib — Cylindric Plane Partition library for rank-3 research.

Re-exports the most commonly used symbols from submodules.
"""

from .partitions import get_part, gen_partitions, gen_compatible, conjugate, partitions_of
from .profiles import profile_to_mu, all_profiles, representative_profiles
from .enumeration import (compute_gk, enumerate_cpps_with_max,
                          extract_z_coefficient, gk_z1_coefficients,
                          compute_g_m, compute_gk_both_stats)
from .qseries import (qpoch, qbinomial, compute_qpoch_m, gk_dict_to_poly,
                       compute_h_m, compute_Q, borodin_product)
from .qsym import (compositions_of, refines, coarsens, q_specialization,
                    M_to_L, M_to_L_bounded)
from .ground_set import frequency_vector, enumerate_all_triples
from .cli import parse_args
