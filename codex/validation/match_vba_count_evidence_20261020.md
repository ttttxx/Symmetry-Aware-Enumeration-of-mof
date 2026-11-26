# VBA 独立层计数原文存档（2026-10-20）

按作者要求，论文中独立层计数与运行限制并回 Results 一句说明。本文件保留原有数学核验细节，供后续续接与审计；不作为论文补充材料。运行记录位于 results/vba_rho_stratification_20261010/。

```latex
\subsection{Independent exact-count and closure checks for \texttt{vba}}
\label{app:vba_checks}

For \texttt{vba}, Lagrange's theorem restricts stabilizer orders to $16,8,4,2,1$. Each run used $B_s=20{,}000$, a global cap of $M_{\mathrm{total}}=100{,}000=\sum_sB_s$, and a 60\,s limit per fixed-point task. No task timed out, and neither run reached the global cap. The quota guarantee therefore classifies each exhausted layer as complete and each remaining layer as quota-filled. The totals 41{,}984 and 60{,}549 reflect separate quotas whose unused slots remain unallocated.
% 中文翻译：对于 \texttt{vba}，拉格朗日定理将稳定子阶数限定为 $16,8,4,2,1$。每次运行均设置 $B_s=20{,}000$、总上限 $M_{\mathrm{total}}=100{,}000=\sum_sB_s$，并为每个固定点任务设置 60 秒时限。没有任务超时，且两次运行均未达到总上限。因此，配额保证可将已耗尽层判为完整，将其余层判为已达到配额。41{,}984 和 60{,}549 的总数来自各层独立配额，未使用的名额不会转配到其他层。

The stratum sizes can also be checked without visiting their colorings. Let $a(H)$ count labeled assignments with stabilizer exactly $H$. Equation~\ref{eq:fixed_layer} partitions $\mathrm{Fix}(H)$ by stabilizer supergroups, yielding the descending recursion
\begin{equation}
\begin{aligned}
a(H)&=\prod_d|\mathcal{C}_d|^{m_{H,d}}-\sum_{K>H}a(K),\\
|\Omega_s|&=\frac{s}{|\Gamma|}\sum_{|H|=s}a(H).
\end{aligned}
\label{supp:eq:exact_layer_count}
\end{equation}
Here the second line divides the number of labeled assignments in the exact-$s$ stratum by their common orbit size $|\Gamma|/s$. For this benchmark action, every nonidentity element is an involution and fixes no vertex. Thus $\Gamma\cong C_2^4$ has two regular vertex orbits, and $m_{H,4}=32/|H|$. Its subgroup counts at orders $16,8,4,2,1$ are $1,15,35,15,1$. Applying Eq.~\ref{supp:eq:exact_layer_count} gives
\begin{equation}
\begin{array}{c|r|r}
s & |\Omega_s|\ (|\mathcal{C}_4|=2) & |\Omega_s|\ (|\mathcal{C}_4|=3)\\ \hline
16 & 4 & 9\\
8 & 90 & 540\\
4 & 1{,}890 & 55{,}440\\
2 & 119{,}880 & 80{,}628{,}480\\
1 & 268{,}375{,}032 & 115{,}813{,}721{,}475{,}072
\end{array}
\label{supp:eq:vba_exact_sizes}
\end{equation}
These integer fixed-set calculations independently reproduce every quota or exhaustion count in Table~\ref{tab:vba_bounded} against the quota guarantee stated in the main text. They certify layer sizes without materializing the corresponding orbit representatives.
% 中文翻译：还可在不逐个访问着色的情况下独立核查层大小。令 $a(H)$ 为稳定子恰好等于 $H$ 的有标签指派数。固定点集按稳定子超群划分，得到公式~\ref{supp:eq:exact_layer_count} 的降序递推；第二行按共同轨道大小将有标签计数换成轨道数。该基准群的非恒等元素均为无不动顶点的对合，因此群同构于 $C_2^4$，顶点分成两个正则轨道，并有 $m_{H,4}=32/|H|$。阶数 $16,8,4,2,1$ 的子群数分别为 $1,15,35,15,1$；递推计算得到公式~\ref{supp:eq:vba_exact_sizes} 的两个候选池完整精确层大小。代入正文配额公式，可独立复现表~\ref{tab:vba_bounded} 的全部配额与耗尽计数，无需生成对应轨道代表列表。


```
