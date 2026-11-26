"""Re-export current MATCH pore descriptors from original Zeo++ Di/Df values."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / 'data_availability/mof-enumeration-reconstruction'
LOCAL_RAW = ROOT / 'results/co2n2_matched_ordering_expansion_2026-10-10/targeted_zeopp'
UPSTREAM_RAW = ROOT.parent / 'CO2_N2/outputs/co2n2_ffc_screening/batch_002_clean_pool_5000/zeopp'


def read_table(path):
    with path.open(newline='') as f:
        reader = csv.DictReader(f)
        return reader.fieldnames, list(reader)


def write_table(path, fields, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def raw_diameters(candidate_id):
    for root in (LOCAL_RAW, UPSTREAM_RAW):
        path = root / candidate_id / (candidate_id + '.res')
        if path.exists():
            di, df, dif = path.read_text().split()[-3:]
            assert float(di) >= float(df)
            return di, df
    raise FileNotFoundError(candidate_id)


def main():
    tables = [
        PACKAGE / 'tables/ffc_co2n2_ordering_check.csv',
        PACKAGE / 'tables/targeted_zeopp_descriptors.csv',
        ROOT / 'results/co2n2_matched_ordering_expansion_2026-10-10/targeted_zeopp_descriptors.csv',
        ROOT / 'results/co2n2_repeat_run_2026-10-19/candidate_metrics_five_run_summary.csv',
    ]
    for path in tables:
        fields, rows = read_table(path)
        pld, lcd = ('pld_angstrom', 'lcd_angstrom') if 'pld_angstrom' in fields else ('pld', 'lcd')
        for row in rows:
            di, df = raw_diameters(row['candidate_id'])
            row[pld], row[lcd] = df, di
        write_table(path, fields, rows)
        print(f'Re-exported {len(rows)} rows: {path.relative_to(ROOT)}')

    _, candidates = read_table(tables[0])
    path = ROOT / 'results/co2n2_repeat_run_2026-10-19/matched_group_ranges_five_run_mean.csv'
    fields, ranges = read_table(path)
    for row in ranges:
        group = [c for c in candidates if c['group_label'] == row['group']]
        assert len(group) == int(row['n'])
        for metric in ('pld_angstrom', 'lcd_angstrom'):
            values = [float(c[metric]) for c in group]
            row[metric + '_min'] = format(min(values), '.10g')
            row[metric + '_max'] = format(max(values), '.10g')
    write_table(path, fields, ranges)
    print('Recomputed nine current group ranges; GCMC fields retained.')


if __name__ == '__main__':
    main()
