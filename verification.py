# mof/verification.py
import math
from functools import reduce
from operator import mul


class BridgeVerifier:
    """
    实现 main.tex 中的 Proposition 2.1: Algorithmic Self-Verification
    Identity: Sum(1/|Stab(s)|) = |Fix(e)| / |G|
    """

    def __init__(self, group_size):
        self.group_size = group_size

    def _nCr_multi(self, n, counts):
        """
        计算多项式系数 (Multinomial Coefficient): n! / (n1! * n2! * ...)
        用于计算固定组分下的总方案数 |Fix(e)|
        """
        denom = reduce(mul, (math.factorial(c) for c in counts), 1)
        return math.factorial(n) // denom

    def compute_rhs_theoretical(self, num_nodes, atom_types, composition=None,
                                 topology_requirements=None, user_selected_ligands=None):
        """
        计算公式右侧 (RHS): Counting Value = |Fix(e)| / |G|
        |Fix(e)| 是所有可能的着色方案总数（即恒等置换的不动点数）。

        当提供 topology_requirements 和 user_selected_ligands 时，
        |Fix(e)| 会受到度约束限制：度数为 d 的节点只能从 user_selected_ligands[d] 中选择配体。
        """
        if composition:
            # Case 1: 固定组分 (Fixed Composition)
            # |Fix(e)| = 全排列数 (Multinomial Coefficient)
            counts = list(composition.values())
            total_schemes = self._nCr_multi(num_nodes, counts)
        elif topology_requirements and user_selected_ligands:
            # Case 2: 带度约束的所有组分
            # |Fix(e)| = Π_{d} (|ligands[d]|)^{count_d}
            total_schemes = 1
            for deg, count in topology_requirements.items():
                n_choices = len(user_selected_ligands.get(deg, atom_types))
                total_schemes *= (n_choices ** count)
        else:
            # Case 3: 无约束的所有组分
            # |Fix(e)| = (颜色数)^节点数
            num_colors = len(atom_types)
            total_schemes = num_colors ** num_nodes

        rhs = total_schemes / self.group_size
        return rhs, total_schemes

    def compute_lhs_constructive(self, results):
        """
        计算公式左侧 (LHS): Constructive Sum = Sum(1 / |Stab(s)|)
        只有当算法找到了所有轨道代表时，这个和才等于 RHS。
        """
        lhs = 0.0
        for res in results:
            stab_size = res.get('stabilizer_size', 0)
            if stab_size > 0:
                lhs += (1.0 / stab_size)
        return lhs

    def verify(self, results, num_nodes, atom_types, composition=None, tolerance=1e-9,
               topology_requirements=None, user_selected_ligands=None):
        """
        执行验证逻辑
        """
        print("\n=== Algorithmic Self-Verification ===")

        # 1. 计算 LHS (构造值)
        lhs = self.compute_lhs_constructive(results)

        # 2. 计算 RHS (理论值)
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