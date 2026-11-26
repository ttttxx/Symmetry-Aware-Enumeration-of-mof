"""Count abstract terms against subsequent English LaTeX source, excluding comments."""
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'paper/main_MATCH.tex'


def english_source(text):
    return re.sub(r'(?<!\\)%[^\n]*', '', text)


def normalize(text):
    text = text.replace('CO$_2$', 'CO2').replace('N$_2$', 'N2')
    text = re.sub(r'\\(?:cite|ref|label|eqref|pageref)\{[^}]*\}', '', text)
    text = re.sub(r'\\(?:begin|end)\{[^}]*\}', ' ', text)
    text = re.sub(r'\\[A-Za-z]+\*?(?:\[[^]]*\])?', '', text)
    text = text.replace('{', '').replace('}', '').replace('~', ' ')
    return re.sub(r'\s+', ' ', text).strip().lower()


TERMS = [
    'metal-organic frameworks', 'MOFs', 'tunable structures', 'gas storage',
    'catalysis', 'separation', 'building-unit arrangements', 'node types',
    'fixed MOF net', 'ordered arrangements', 'inequivalent', 'selected symmetries',
    'periodic models', 'compatible building-unit connections and orientations',
    'degree-compatible node assignments', 'quotient-graph colorings',
    'fixed linker template', 'orientation-preserving', 'candidate equivalence',
    'subgroup-lattice fixed-point searches', 'subgroup fixed-point searches',
    'canonicalization', 'orbit representatives', 'without duplication',
    'Kabsch-based node alignment', 'node alignment', 'linker placement',
    'oriented periodic models', 'selected representatives',
    'geometric construction audit', 'target topology', 'short-contact check',
    'reloadable', 'crystallographic information files', 'CIFs', 'CO2/N2',
    'adsorption screening', 'single-point density functional theory', 'DFT',
    'energy comparisons', 'composition groups',
]


def count(text, term):
    return len(re.findall(r'(?<![a-z0-9])' + re.escape(term.lower()) + r'(?![a-z0-9])', text))


def main():
    raw = SOURCE.read_text()
    clean = english_source(raw)
    abstract = clean.split(r'\begin{abstract}', 1)[1].split(r'\end{abstract}', 1)[0]
    later = clean.split(r'\end{abstract}', 1)[1]
    main_text, si = later.split(r'\appendix', 1)
    zones = {'abstract': normalize(abstract), 'main_text': normalize(main_text), 'SI': normalize(si)}
    rows = []
    for term in TERMS:
        row = {'term': term}
        row.update({key: count(text, term) for key, text in zones.items()})
        row['later_total'] = row['main_text'] + row['SI']
        row['whole_total'] = row['abstract'] + row['later_total']
        rows.append(row)
    out = ROOT / 'codex/validation/match_abstract_term_counts_20261111.csv'
    with out.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    words = re.findall(r'[a-z]+(?:-[a-z]+)*', zones['abstract'])
    abstract_counts = Counter(words)
    later_counts = Counter(re.findall(r'[a-z]+(?:-[a-z]+)*', zones['main_text'] + ' ' + zones['SI']))
    one_off = [w for w, n in abstract_counts.items() if n == 1 and later_counts[w] == 0]
    details = {
        'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'scope': 'English source after abstract, including prose, headings, captions and algorithms; main text and SI separated; comments excluded; rendered bibliography not included.',
        'matching': 'Case insensitive exact word/phrase; singular/plural and synonyms NOT combined; LaTeX commands removed; CO2/N2 normalized.',
        'abstract_tokens': len(words),
        'abstract_word_forms_absent_later': one_off,
        'counts': rows,
    }
    (out.with_suffix('.json')).write_text(json.dumps(details, indent=2) + '\n')
    for row in rows:
        print(f"{row['term']}: abstract {row['abstract']}, main {row['main_text']}, SI {row['SI']}")
    print('Word forms occurring once in abstract and nowhere later:', ', '.join(one_off))


if __name__ == '__main__':
    main()
