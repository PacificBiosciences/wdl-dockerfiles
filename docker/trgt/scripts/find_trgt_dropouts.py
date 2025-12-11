#!/usr/bin/env python3
"""
Identify TRGT loci that are affected by haplotype-specific coverage dropouts.

This script parses a TRGT spanning reads BAM file with haplotype tags (HP) and read quality tags (rq),
and a BED4+ TRGT repeat catalog of TRGT regions of interest. It counts the number of reads per haplotype
spanning each TRGT region, and determines whether there are coverage dropouts based on user-defined
thresholds, outputting a tab-delimited summary for tandem repeat regions with detected dropouts.
"""

__version__ = '0.2.0'

import argparse
import gzip
import logging
import math
import sys
from typing import Dict, Optional

import pysam

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y%m%dT%H:%M:%S%z', level=logging.INFO)


def parse_sex(sex_str: Optional[str]) -> str:
  """Parse"""
  if sex_str is None:
    logging.warning('Sex not provided; assuming FEMALE for dropout determination.')
    return 'FEMALE'
  if sex_str not in ('MALE', 'FEMALE'):
    logging.warning(f'Unrecognized sex: {str(sex_str)}; assuming FEMALE for dropout determination.')
    return 'FEMALE'
  return sex_str


def parse_tr_id(label: str) -> str:
  """Extract TR id from a label string containing 'ID=...'. Raises ValueError if not found."""
  for part in label.split(';'):
    if part.startswith('ID='):
      return part.split('=', 1)[1]
  raise ValueError(f'Missing ID=... in label: {label}')


class TandemRepeatRegion:
  """Represent a tandem repeat region with per-haplotype read counts."""

  def __init__(self, chrom: str, chrom_start: int, chrom_end: int, label: str) -> None:
    self.chrom: str = chrom
    self.chrom_start: int = int(chrom_start)
    self.chrom_end: int = int(chrom_end)
    self.label: str = label
    self.trid: str = parse_tr_id(label)
    self.reads: Dict[str, int] = {'0': 0, '1': 0, '2': 0, 'fail': 0}

  def add_haplotagged_read(self, key) -> None:
    """Add a read to the appropriate haplotype bucket."""
    self.reads[key] += 1


def rq_from_bq(bq: 'list[int]') -> float:
  """Compute read quality from an array of phred scaled base qualities"""
  read_len = len(bq)
  expectedErrors = sum([math.pow(10, -0.1 * x) for x in bq])
  return expectedErrors / read_len if read_len > 0 else 0.0


def load_bed(trbed: str, chrom: Optional[str] = None) -> Dict[str, TandemRepeatRegion]:
  """Load BED4+ regions and return dict mapping TR ID -> TandemRepeatRegion."""

  def open_bed(path: str):
    if path.endswith('.gz'):
      return gzip.open(path, 'rt')
    return open(path, 'r')

  tandem_repeats: Dict[str, TandemRepeatRegion] = {}
  with open_bed(trbed) as fh:
    for line in fh:
      if not line.strip() or line.startswith('#'):
        continue
      cols = line.rstrip('\n').split('\t')
      if len(cols) < 4:
        logging.warning(f'Skipping malformed BED line: {line.rstrip()}')
        continue
      if chrom and cols[0] != chrom:
        continue
      tr = TandemRepeatRegion(cols[0], int(cols[1]), int(cols[2]), cols[3])
      tandem_repeats[tr.trid] = tr
  logging.info(
    f'Loaded {len(tandem_repeats)} tandem repeat regions from {trbed}' + (f' for chrom {chrom}' if chrom else '')
  )
  return tandem_repeats


def assign_haplotags(
  tandem_repeats: Dict[str, TandemRepeatRegion], spanning_bam: str, chrom: Optional[str] = None
) -> None:
  """Walk haplotagged spanning BAM and assign read counts to TR regions by haplotype."""
  with pysam.AlignmentFile(spanning_bam, 'rb') as tb:
    for r in tb.fetch(chrom) if chrom else tb.fetch():
      if r.has_tag('TR') and r.get_tag('TR') in tandem_repeats:
        if r.has_tag('rq'):
          rq = float(r.get_tag('rq'))
        else:
          logging.warning(
            f'Read {str(r.query_name) if r.query_name else ""} missing rq tag; computing from base qualities.'
          )
          rq = rq_from_bq(list(r.query_qualities)) if r.query_qualities else 0.0
        if rq < 0.99:
          hp = 'fail'
        else:
          hp = str(r.get_tag('HP')) if r.has_tag('HP') else '0'
        tandem_repeats[str(r.get_tag('TR'))].add_haplotagged_read(hp)


def determine_dropout(
  hap1: int, hap2: int, unphased: int, fail: int, chrom: str, sex: str, coverage: int, trid: str
) -> str:
  """Determine dropout label given counts, chrom, sex and coverage threshold."""
  if sex == 'FEMALE' and chrom in ('chrY', 'Y'):
    # no expectation of coverage on chrY for females
    return ''
  if hap1 >= coverage and hap2 >= coverage:
    # both haplotypes have adequate coverage -> no dropout
    return ''
  if (hap1 < coverage and hap2 >= coverage) or (hap2 < coverage and hap1 >= coverage):
    # one haplotype is low and the other is adequate, call haplotype dropout
    # males are hemizygous for chrX/Y expectation -> do not call haplotype dropout
    if sex == 'MALE' and chrom in ('chrX', 'X', 'chrY', 'Y'):
      return ''
    return 'HaplotypeDropout'
  if hap1 + hap2 + unphased + fail < 2 * coverage:
    return 'FullDropout'
  if hap1 < coverage and hap2 < coverage:
    # both haplotypes are low but unphased is adequate, call phasing dropout
    # males are hemizygous for chrX/Y expectation -> do not call phasing dropout
    if sex == 'MALE' and chrom in ('chrX', 'X', 'chrY', 'Y'):
      return ''
    return 'PhasingDropout'
  raise ValueError(
    f'Unhandled dropout determination case: trid={trid}, hap1={hap1}, hap2={hap2}, unphased={unphased}, fail={fail}, chrom={chrom}, sex={sex}, coverage={coverage}'
  )


def compute_and_print_results(
  tandem_repeats: Dict[str, TandemRepeatRegion], sex: str, coverage: int, print_all: bool
) -> None:
  """Compute dropouts and write tab-delimited output to stdout."""
  print(
    '\t'.join(
      [
        'chrom',
        'start',
        'end',
        'trid',
        'hap1_count',
        'hap2_count',
        'unphased_count',
        'fail_read_count',
        'dropout',
      ]
    ),
    file=sys.stdout,
  )
  for trid, tr in tandem_repeats.items():
    unphased = tr.reads['0']
    hap1 = tr.reads['1']
    hap2 = tr.reads['2']
    fail = tr.reads['fail']
    dropout = determine_dropout(hap1, hap2, unphased, fail, tr.chrom, sex, coverage, trid)
    if not print_all and dropout == '' and fail == 0:
      # no dropout and no failed reads, skip output
      continue
    print(
      f'{tr.chrom}\t{tr.chrom_start}\t{tr.chrom_end}\t{tr.label}\t{hap1}\t{hap2}\t{unphased}\t{fail}\t{dropout}',
      file=sys.stdout,
    )


def main(argv=None):
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('trbed', help='Regions of interest, BED4+ (chrom, chromStart, chromEnd, label)')
  parser.add_argument('spanning_bam', help='Spanning reads bam.')
  parser.add_argument('--chrom', help='Chromosome to analyze', default=None)
  parser.add_argument('--sex', choices=['MALE', 'FEMALE'], default=None)
  parser.add_argument('--coverage', type=int, default=2, help='Minimum per-haplotype coverage for region.')
  parser.add_argument('--print-all', action='store_true', help='Print all regions, not just those with dropouts.')
  parser.add_argument('--version', action='version', version=f'%(prog)s (version {__version__})')
  args = parser.parse_args(argv)

  logging.info(f'Starting find_trgt_dropouts (version {__version__})')

  tandem_repeats = load_bed(args.trbed, args.chrom)
  assign_haplotags(tandem_repeats, args.spanning_bam, args.chrom)
  compute_and_print_results(tandem_repeats, parse_sex(args.sex), args.coverage, args.print_all)
  logging.info('Completed processing')


if __name__ == '__main__':
  main()
