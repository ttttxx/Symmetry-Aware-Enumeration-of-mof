"""Apply approved fresh benchmark data and update the local reproducibility package."""
from pathlib import Path
import collections
import csv
import hashlib
import json
import shutil
import statistics

ROOT=Path(__file__).resolve().parents[1]
PAPER=ROOT/'paper/main_MATCH.tex'
REFRESH=ROOT/'validation/MATCH_UNIFIED_REFRESH_2026-11-23.json'
VBA=ROOT/'validation/MATCH_UNIFIED_VBA_2026-11-23.json'
MANIFEST=ROOT/'validation/MATCH_UNIFIED_MANIFEST_2026-11-23.json'
PACKAGE=ROOT/'data_availability/mof-enumeration-reconstruction'


def number(n):
    return f'{n:,}'


def time_s(t):
    return f'{t:.4g}'


def replace_table(text,label,new):
    position=text.index('\\label{'+label+'}')
    start=text.rfind('\\begin{table}',0,position)
    end=text.index('\\end{table}',position)+len('\\end{table}')
    return text[:start]+new+text[end:]


def table(caption,chinese,label,columns,header,rows):
    return ('\\begin{table}[H]\n\\centering\\scriptsize\n'
        '\\setlength{\\tabcolsep}{1.2mm}\n'
        '\\caption{'+caption+'}\n% 中文翻译：'+chinese+'\n'
        +(('\\label{'+label+'}\n') if label else '')
        +'\\begin{tabular}{'+columns+'}\n\\toprule\n'+header+' \\\\\n\\midrule\n'
        +'\n'.join(' & '.join(row)+' \\\\' for row in rows)
        +'\n\\bottomrule\n\\end{tabular}\n\\end{table}')


if __name__=='__main__':
    data=json.loads(REFRESH.read_text());vba=json.loads(VBA.read_text());manifest=json.loads(MANIFEST.read_text())
    for row in manifest['files']:
        assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    assert len(data['records'])==108 and len(data['summary'])==9
    assert all(row['status']=='exhaustive' for row in data['records'])
    assert all(row['matched_record_digests'] for row in data['summary'])
    assert len(vba['records'])==12 and len(vba['summary'])==2
    summary=data['summary']
    rows=[]
    for row in summary:
        rows.append(['\\texttt{'+row['topology']+'}',f"{row['c3']}/{row['c4']}",number(row['naive']),
            number(row['orbits']),time_s(row['ours_median_s']),time_s(row['casm_records_median_s'])])
    text=PAPER.read_text()
    main_caption=('\\textbf{Matched-group enumeration on benchmark quotient graphs.} Both workflows use the same '
        '$\\Gamma$ and return complete assignments, stabilizer orders, and orbit sizes. '
        'Times are enumeration wall-time medians of three single-threaded runs. '
        'Input pools and timing boundaries are specified in Supplementary Section~\\ref{app:enum_benchmark}.')
    text=replace_table(text,'tab:casm',table(main_caption,
        '基准商图上的同群枚举。两种流程使用相同的$\\Gamma$，并返回完整指派、稳定子阶数及轨道大小。'
        '时间为三次单线程运行的枚举壁钟中位数；输入池和计时边界见补充材料第~\\ref{app:enum_benchmark}节。',
        'tab:casm','llrrrr','Net & $|\\mathcal C_3|/|\\mathcal C_4|$ & $|X|$ & Orbits & This work (s) & CASM (s)',rows[:8]))
    extra_caption='Matched-group enumeration with complete representative records. Times are enumeration wall times; CASM uses the same $\\Gamma$.'
    text=replace_table(text,'tab:matched_extra',table(extra_caption,
        '保留完整代表记录的同群枚举。时间为枚举壁钟耗时，CASM采用相同的$\\Gamma$。',
        'tab:matched_extra','llrrrr','Net & $|\\mathcal C_3|/|\\mathcal C_4|$ & $|X|$ & Orbits & This work (s) & CASM (s)',rows))
    # Keep the existing interpretation pending item 3, and repair its table reference.
    old='Its native full factor group $F$ includes orientation-reversing operations in addition to the orientation-preserving operations retained here.'
    text=text.replace(old,old[:-1]+' (Supplementary Table~\\ref{tab:orbit_group_counts}).')
    text=text.replace('其原生完整因子群 $F$ 还包含本文保向操作之外的反向操作。',
        '其原生完整因子群 $F$ 还包含本文保向操作之外的反向操作（补充表~\\ref{tab:orbit_group_counts}）。')
    text=text.replace('The additional enumeration benchmarks used the ABI-3 implementation snapshot,',
        'The enumeration benchmarks used the current ABI-4 implementation snapshot,')
    text=text.replace('新增枚举基准采用ABI-3实现快照','枚举基准采用当前ABI-4实现快照')
    text=text.replace('with cardinalities given in Table~\\ref{tab:matched_extra}.',
        'with cardinalities given in Table~\\ref{tab:casm} and Supplementary Table~\\ref{tab:matched_extra}.')
    text=text.replace('池大小见表~\\ref{tab:matched_extra}',
        '池大小见表~\\ref{tab:casm}及补充表~\\ref{tab:matched_extra}')
    # The bounded runs use prefixes of the same lexical master list, not historical pools.
    needle='Table~\\ref{tab:vba_bounded} compares per-stratum quotas with DFS prefixes containing the same total number of representatives.'
    text=text.replace(needle,needle+' The ordered label pools and run settings are specified in Supplementary Section~\\ref{app:enum_benchmark}.')
    text=text.replace('表~\\ref{tab:vba_bounded} 将分层配额与代表总数相同的 DFS 前缀进行比较。',
        '表~\\ref{tab:vba_bounded} 将分层配额与代表总数相同的 DFS 前缀进行比较。有序标签池与运行设置见补充材料第~\\ref{app:enum_benchmark}节。')
    needle='These benchmarks enumerated assignments without constructing atomic models.'
    protocol=('The bounded \\texttt{vba} runs used the first two or three labels of the same ordered $\\mathcal C_4$ master list. '
        'Canonical comparisons used label-string lexicographic order, and DFS visited labels in the listed order. '
        'Both policies used the ABI-4 C++ backend with the representative index and completed-prefix cache enabled. '
        'Quota runs used $B_s=20{,}000$, a total cap of 100,000, and a 60\\,s task limit; '
        'DFS prefixes stopped at the corresponding quota run\'s total output. '
        'No task timed out. Layer times are medians of three independent runs and include the work performed within each layer; '
        'initial subgroup setup is outside those layer intervals.\n'
        '% 中文翻译：受限vba运行使用同一有序四连接主列表的前两或前三个标签。规范比较按标签字符串字典序，DFS按表列顺序遍历。'
        '两种策略均使用ABI-4 C++后端并开启代表索引与已完成前缀缓存。配额运行采用每层20,000、总上限100,000及每任务60秒；'
        'DFS前缀在对应配额运行的总输出数处停止。没有任务超时。层耗时为三次独立运行中位数，包含层内执行的工作；初始子群准备在层时间之外。\n\n')
    assert text.count(needle)==1
    text=text.replace(needle,protocol+needle)
    # Group-count table: values are freshly rechecked even though cardinalities are unchanged.
    count_rows=[]
    for row in summary:
        full=next(r['orbits'] for r in data['records'] if r['topology']==row['topology'] and r['c3']==row['c3'] and r['c4']==row['c4'] and r['side']=='casm_full')
        count_rows.append(['\\texttt{'+row['topology']+'}',f"{row['c3']}/{row['c4']}",number(row['orbits']),number(full)])
    old_caption='\\caption{Assignment-orbit counts under the orientation-preserving group $\\Gamma$ and the native full factor group $F$.}'
    position=text.index(old_caption)
    start=text.rfind('\\begin{table}',0,position);end=text.index('\\end{table}',position)+len('\\end{table}')
    text=text[:start]+table('Assignment-orbit counts under the orientation-preserving group $\\Gamma$ and the native full factor group $F$.',
        '保向群$\\Gamma$与原生完整因子群$F$下的指派轨道计数。','tab:orbit_group_counts','llrr',
        'Net & $|\\mathcal C_3|/|\\mathcal C_4|$ & Orbits ($\\Gamma$) & Orbits ($F$)',count_rows)+text[end:]
    range_rows=[];core_rows=[]
    for row in summary:
        runs=[r for r in data['records'] if (r['topology'],r['c3'],r['c4'])==(row['topology'],row['c3'],row['c4'])]
        ours=[r for r in runs if r['side']=='ours'];casm=[r for r in runs if r['side']=='casm_records']
        prefix=['\\texttt{'+row['topology']+'}',f"{row['c3']}/{row['c4']}"]
        range_rows.append(prefix+[time_s(min(r['time_s'] for r in ours))+'--'+time_s(max(r['time_s'] for r in ours)),
            time_s(min(r['time_s'] for r in casm))+'--'+time_s(max(r['time_s'] for r in casm)),
            f"{row['ours_median_rss_mib']:.1f}",f"{row['casm_records_median_rss_mib']:.1f}"])
        core_rows.append(prefix+[time_s(statistics.median(r['native_kernel_s'] for r in ours)),time_s(row['casm_count_median_s'])])
    for caption,zh,cols,header,rr in [
        ('Observed wall-time ranges over three runs and peak process memory medians. Time ranges are in seconds and memory is in MiB.',
         '三次运行的壁钟耗时范围及进程峰值内存中位数。耗时单位为秒，内存单位为MiB。','llrrrr',
         'Net & Pool & Time range (ours) & Time range (CASM) & RSS (ours) & RSS (CASM)',range_rows),
        ('Native-core time of this work and a separate CASM count-only comparison under $\\Gamma$ (seconds, medians of three runs). Native-core time is a subset of enumeration wall time.',
         '本文原生核时间与$\\Gamma$下CASM仅计数的独立对照，单位为秒，取三次中位数。原生核时间是枚举壁钟时间的一部分。','llrr',
         'Net & Pool & Native core (ours) & CASM count only',core_rows)]:
        position=text.index('\\caption{'+caption+'}')
        start=text.rfind('\\begin{table}',0,position);end=text.index('\\end{table}',position)+len('\\end{table}')
        text=text[:start]+table(caption,zh,None,cols,header,rr)+text[end:]
    largest=next(r for r in summary if (r['topology'],r['c3'],r['c4'])==('fjh',1,4))
    text=text.replace('887.5\\,MiB for this work and 720.1\\,MiB for CASM',
        f"{largest['ours_median_rss_mib']:.1f}\\,MiB for this work and {largest['casm_records_median_rss_mib']:.1f}\\,MiB for CASM")
    text=text.replace('887.5 MiB和720.1 MiB',f"{largest['ours_median_rss_mib']:.1f} MiB和{largest['casm_records_median_rss_mib']:.1f} MiB")
    lil=next(r for r in summary if r['topology']=='lil')
    text=text.replace('1.952\\,s for this work and 6.817\\,s',f"{time_s(lil['ours_median_s'])}\\,s for this work and {time_s(lil['casm_records_median_s'])}\\,s")
    text=text.replace('分别为1.952秒和6.817秒',f"分别为{time_s(lil['ours_median_s'])}秒和{time_s(lil['casm_records_median_s'])}秒")
    # Keep every bounded count and composition-vector value tied to the new runs.
    position=text.index('\\label{tab:vba_bounded}')
    start=text.index('\\begin{tabular}',position);end=text.index('\\end{tabular}',start)
    body='\\begin{tabular}{llccc}\n\\toprule\n\\textbf{Pool} & $\\boldsymbol{\\rho}$ & \\textbf{DFS prefix} & \\textbf{Per-$\\rho$ quota} & \\textbf{Layer time (s)} \\\\\n\\midrule\n'
    for item in vba['summary']:
        c4=item['c4'];body+='\\multicolumn{5}{l}{$|\\mathcal{C}_4|='+str(c4)+'$; matched total = '+number(item['quota']['returned'])+'} \\\\\n'
        for order in [16,8,4,2,1]:
            pre=item['prefix']['distribution'][str(order)];quota=item['quota']['distribution'][str(order)]
            t=next(r['median_s'] for r in item['quota']['layers'] if r['stabilizer_order']==order)
            body+=' & '.join([f'$|\\mathcal C_4|={c4}$',f'${order/16:g}$',
                f"{number(pre['returned'])}\\;({pre['composition_vectors']})",f"{number(quota['returned'])}\\;({quota['composition_vectors']})",f'{t:.4f}'])+' \\\\\n'
        if c4==2:body+='\\addlinespace\n'
    text=text[:start]+body+'\\bottomrule\n'+text[end:]
    text=text.replace("layer time is the time spent on the quota schedule's fixed-point tasks.",
        'layer times are medians of three runs for the layer intervals defined in Supplementary Section~\\ref{app:enum_benchmark}.')
    text=text.replace('层耗时为配额调度在该层固定点任务上的耗时。',
        '层耗时为三次运行中位数，其层区间定义见补充材料第~\\ref{app:enum_benchmark}节。')
    old='source tables for the CASM comparison and DFT energies, CO$_2$/N$_2$ screening data, bounded-\\texttt{vba} records, and the X-site fusion audit.'
    text=text.replace(old,'source tables for the matched-group CASM benchmarks and DFT energies, CO$_2$/N$_2$ screening data, bounded-\\texttt{vba} records, and the X-site fusion audit. Versioned benchmark records, input pools, and source snapshots document the reported enumeration results.')
    text=text.replace('CASM 对比和 DFT 能量源数据表、CO$_2$/N$_2$ 筛选数据、受限 \\texttt{vba} 记录及 X 位点融合审计数据。',
        '同群CASM基准和DFT能量源数据表、CO$_2$/N$_2$筛选数据、受限\\texttt{vba}记录及X位点融合审计数据。版本化基准记录、输入标签池和源码快照记录了本文枚举结果。')
    PAPER.write_text(text)
    # Freeze final manuscript separately: no executable source changes after measurement.
    snap=ROOT/manifest['snapshot_root'];shutil.copy2(PAPER,snap/'main_MATCH_after.tex')
    target=PACKAGE/'reproducibility/unified_20261123'
    shutil.copytree(snap,target/'snapshot',dirs_exist_ok=True)
    for source in [REFRESH,REFRESH.with_suffix('.jsonl'),REFRESH.with_suffix('.csv'),VBA,VBA.with_suffix('.jsonl'),MANIFEST,Path(__file__)]:
        dest=target/source.name;shutil.copy2(source,dest)
    (target/'README.md').write_text('''# Current enumeration benchmarks (2026-11-23)

The manuscript enumeration tables use this ABI-4 C++/Trie snapshot and the lexical template pools in the measurement JSON. All 108 exhaustive measurements completed and all normalized complete records matched CASM under the same Gamma. Twelve bounded vba runs use prefixes of the same C4 master list. Quota and prefix runs have equal output totals; each policy was repeated three times in independent processes.

Reproduce from the project root represented by snapshot/: build the native library with tools/build_enumeration_core.py, then run experiment/match_refresh_compare.py --skip-legacy --output NEW.json using a Python environment with psutil, and validation/benchmark_unified_vba_20261123.py --output NEW_VBA.json with Python 3.11. The matched harness defaults to python3.11 workers. Dependencies include NumPy, SymPy, pymatgen, libcasm-configuration 2.3.1 and psutil. The bundled dylib is macOS arm64; rebuild on another platform.

Enumeration wall time excludes imports, graph/group construction and verification; it includes subgroup tasks, search, record storage, decoding and sorting. Native-core time is a subset. Layer intervals exclude the initial subgroup setup. Old dated CSV files elsewhere in the package are historical provenance and are not the source of current manuscript enumeration timings. This local package has not been uploaded or pushed by this task.
''')
    mapping=PACKAGE/'figure_table_mapping.md';m=mapping.read_text()
    m=m.replace('`tables/vba_rho_stratification_20261010.csv` |',
        '`reproducibility/unified_20261123/MATCH_UNIFIED_VBA_2026-11-23.json` (per-layer medians, counts, vectors and ordered digests) |')
    m=m.replace('`tables/casm_enumeration_comparison.csv` |',
        '`reproducibility/unified_20261123/MATCH_UNIFIED_REFRESH_2026-11-23.json` and `.csv` (current ABI-4, matched Gamma and records) |')
    m+='\n## Current supplementary tables\n\n| Manuscript item | Package files |\n|---|---|\n| Matched-group records, Gamma/F counts, time ranges/RSS, native core and CASM count-only tables | `reproducibility/unified_20261123/MATCH_UNIFIED_REFRESH_2026-11-23.json` |\n| G1–G9 node compositions | `tables/ffc_co2n2_group_compositions.csv`, checked against all 30 `ffc_co2n2_ordering_check.csv` records |\n| Enumeration implementation and exact input version | `reproducibility/unified_20261123/MATCH_UNIFIED_MANIFEST_2026-11-23.json`, `snapshot/` |\n'
    mapping.write_text(m)
    readme=PACKAGE/'README.md';r=readme.read_text().replace('From Topological Design to Pre-Relaxation MOF Candidate Models: Symmetry-Aware Enumeration and Reconstruction',
        'Symmetry-Constrained Enumeration and Oriented Construction of Ordered MOF Candidates')
    r+='\n## Current enumeration source of truth (2026-11-23)\n\nThe current main and supplementary enumeration tables use `reproducibility/unified_20261123/`, with the ABI-4 implementation and one common lexicographic template-pool definition. The older CASM and vba CSV files listed above are dated provenance; their timings are not used in the current manuscript. The 72-model audit, adsorption and DFT samples retain their separately documented construction pools. Group compositions for the adsorption table are in `tables/ffc_co2n2_group_compositions.csv`. No remote upload was performed.\n'
    readme.write_text(r)
    versions=PACKAGE/'code_version.md';versions.write_text(versions.read_text()+'\n## Current ABI-4 enumeration replacement (2026-11-23)\n\nAll manuscript enumeration timings now come from the current-source rerun in `reproducibility/unified_20261123/`. Exact source/library/input hashes are in its manifest, with the measured files preserved in `snapshot/`. Both matched-group and bounded-output experiments use the same lexical template-pool convention. The default CLI/API backend is now auto (C++ preferred, optimized Python fallback); benchmarks explicitly required C++. The new local files have not been pushed or verified as public by this task.\n')
    manifest_file=PACKAGE/'manifest.csv'
    previous={r['package_path']:r for r in csv.DictReader(manifest_file.open())}
    entries=[]
    for path in sorted(PACKAGE.rglob('*')):
        if not path.is_file() or path==manifest_file or '__pycache__' in path.parts:continue
        relative=str(path.relative_to(PACKAGE))
        if relative in previous:entries.append(previous[relative])
        else:entries.append(dict(package_path=relative,source_path='Local verified current-source preparation 2026-11-23',description='Current benchmark measurement, frozen source/input or group-composition mapping; see figure_table_mapping.md'))
    with manifest_file.open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=['package_path','source_path','description']);writer.writeheader();writer.writerows(entries)
    print('Updated current main/supplementary times and counts; packaged',len(entries),'files')
