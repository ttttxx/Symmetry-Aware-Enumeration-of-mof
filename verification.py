import math
from functools import reduce
from operator import mul

class BridgeVerifier:
    """Check the constructive sum against the corresponding counting identity."""

    def __init__(self, group_size):
        self.group_size = group_size

    def _nCr_multi(self, n, counts):
        """Return the multinomial coefficient for a fixed composition."""
        denom = reduce(mul, (math.factorial(c) for c in counts), 1)
        return math.factorial(n) // denom

    def compute_rhs_theoretical(self, num_nodes, atom_types, composition=None,
                                 topology_requirements=None, user_selected_ligands=None):
        if composition:
            counts = list(composition.values())
            total_schemes = self._nCr_multi(num_nodes, counts)
        elif topology_requirements and user_selected_ligands:
            total_schemes = 1
            for deg, count in topology_requirements.items():
                n_choices = len(user_selected_ligands.get(deg, atom_types))
                total_schemes *= (n_choices ** count)
        else:
            num_colors = len(atom_types)
            total_schemes = num_colors ** num_nodes

        rhs = total_schemes / self.group_size
        return rhs, total_schemes

    def compute_lhs_constructive(self, results):
        lhs = 0.0
        for res in results:
            stab_size = res.get('stabilizer_size', 0)
            if stab_size > 0:
                lhs += (1.0 / stab_size)
        return lhs

    def verify(self, results, num_nodes, atom_types, composition=None, tolerance=1e-9,
               topology_requirements=None, user_selected_ligands=None):
        print("\n=== Algorithmic Self-Verification ===")

        lhs = self.compute_lhs_constructive(results)

        rhs, total_schemes = self.compute_rhs_theoretical(
            num_nodes, atom_types, composition,
            topology_requirements=topology_requirements,
            user_selected_ligands=user_selected_ligands)

        diff = abs(lhs - rhs)
        passed = diff < tolerance

        print(f"  Total Schemes |Fix(e)|: {total_schemes}")
        print(f"  Group Size    |G|     : {self.group_size}")
        print(f"  ----------------------------------------")
        print(f"  LHS (Constructive Sum): {lhs:.10f}")
        print(f"  RHS (Counting Value)  : {rhs:.10f}")
        print(f"  Diff                  : {diff:.10e}")

        if passed:
            print("  [SUCCESS] Identity Verified! The enumeration is Complete.")
        else:
            print("  [FAILURE] Mismatch detected. The enumeration is likely incomplete.")

        return passed, {
            "LHS": lhs,
            "RHS": rhs,
            "diff": diff,
            "verified": passed
        }
