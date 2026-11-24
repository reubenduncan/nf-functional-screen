#!/usr/bin/env python3
"""Attach DIAMOND evidence to coordinates in Prodigal protein headers."""
import argparse
import csv
import json
import re
from collections import Counter

p = argparse.ArgumentParser()
for arg in ('sample', 'proteins', 'ncyc', 'bacmet', 'ncyc-map', 'output', 'summary'):
    p.add_argument('--' + arg, required=True)
a = p.parse_args()
with open(a.ncyc_map) as handle:
    mapping = {}
    for line in handle:
        fields = line.rstrip('\n').split('\t')
        if len(fields) >= 2:
            mapping.setdefault(fields[0], set()).add(fields[1])

columns = ['sample', 'contig', 'start_bp', 'end_bp', 'strand', 'gene_id',
           'prodigal_attributes', 'database', 'evidence_class', 'gene_label',
           'reference_id', 'identity_pct', 'alignment_aa', 'query_length_aa',
           'reference_length_aa', 'query_coverage_pct', 'reference_coverage_pct',
           'query_start_aa', 'query_end_aa', 'reference_start_aa', 'reference_end_aa',
           'evalue', 'bitscore', 'reference_description']
by_query = {}
for source, filename in [('ncyc', a.ncyc), ('bacmet', a.bacmet)]:
    with open(filename) as handle:
        for fields in csv.reader(handle, delimiter='\t'):
            if len(fields) != 13:
                raise ValueError(f'{filename}: expected 13 DIAMOND columns, got {len(fields)}')
            q, ref, ident, alen, qlen, slen, qs, qe, ss, se, ev, bits, desc = fields
            if source == 'ncyc':
                label = ';'.join(sorted(mapping.get(ref, set())))
                evidence = 'nitrogen_cycle_candidate' if label else 'non_target_homologue'
            else:
                label = ref.split('|')[1] if '|' in ref else ref
                evidence = 'experimentally_supported_reference_homologue'
            record = dict(sample=a.sample, gene_id=q, database=source, evidence_class=evidence,
                          gene_label=label, reference_id=ref, identity_pct=ident, alignment_aa=alen,
                          query_length_aa=qlen, reference_length_aa=slen,
                          query_coverage_pct=f'{100 * (int(qe)-int(qs)+1)/int(qlen):.2f}',
                          reference_coverage_pct=f'{100 * (int(se)-int(ss)+1)/int(slen):.2f}',
                          query_start_aa=qs, query_end_aa=qe, reference_start_aa=ss,
                          reference_end_aa=se, evalue=ev, bitscore=bits, reference_description=desc)
            by_query.setdefault(q, []).append(record)

found = set()
counts = Counter()
with open(a.output, 'w', newline='') as out:
    writer = csv.DictWriter(out, columns, delimiter='\t')
    writer.writeheader()
    with open(a.proteins) as handle:
        for line in handle:
            if not line.startswith('>'):
                continue
            parts = [x.strip() for x in line[1:].rstrip('\n').split('#')]
            gene = parts[0].split()[0]
            if gene not in by_query:
                continue
            match = re.fullmatch(r'(.+)_([0-9]+)', gene)
            if not match or len(parts) < 5 or parts[3] not in ('1', '-1'):
                raise ValueError('Expected Prodigal protein header with coordinates: ' + line.strip())
            start, end = int(parts[1]), int(parts[2])
            if start < 1 or end < start:
                raise ValueError('Invalid Prodigal coordinates: ' + line.strip())
            for record in by_query[gene]:
                record.update(contig=match.group(1), start_bp=start, end_bp=end,
                              strand='+' if parts[3] == '1' else '-', prodigal_attributes=parts[4])
                writer.writerow(record)
                counts[record['database'] + ':' + record['evidence_class']] += 1
            found.add(gene)
missing = set(by_query) - found
if missing:
    raise ValueError(f'{len(missing)} aligned protein identifiers absent from FASTA; example: {next(iter(missing))}')
with open(a.summary, 'w') as out:
    json.dump({'sample': a.sample, 'aligned_proteins': len(found), 'hits_by_database_and_class': dict(counts)}, out, indent=2)
    out.write('\n')
