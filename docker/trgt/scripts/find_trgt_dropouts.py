#!/usr/bin/env python3
"""
Identify TRGT loci that are affected by haplotype-specific coverage dropouts.

This script parses a TRGT spanning reads BAM file with haplotype tags (HP) and read quality tags (rq),
and a BED4+ TRGT repeat catalog of TRGT regions of interest. It counts the number of reads per haplotype
spanning each TRGT region, and determines whether there are coverage dropouts based on user-defined
thresholds, outputting a tab-delimited summary for tandem repeat regions with detected dropouts.
"""

__version__ = '0.3.0'

import argparse
import gzip
import logging
import math
import sys
from typing import Dict, Tuple, Optional

import pysam

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y%m%dT%H:%M:%S%z', level=logging.INFO)


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

  def __repr__(self) -> str:
    return f'TandemRepeatRegion(chrom={self.chrom}, start={self.chrom_start}, end={self.chrom_end}, trid={self.trid}, reads={self.reads})'

  def add_haplotagged_read(self, key) -> None:
    """Add a read to the appropriate haplotype bucket."""
    self.reads[key] += 1


def rq_from_bq(bq: 'list[int]') -> float:
  """Compute read quality from an array of phred scaled base qualities"""
  read_len = len(bq)
  expectedErrors = sum([math.pow(10, -0.1 * x) for x in bq])
  return 1 - (expectedErrors / read_len) if read_len > 0 else 0.0


def open_bed(path: str):
  if path.endswith('.gz'):
    return gzip.open(path, 'rt')
  return open(path, 'r')


def load_trgt_catalog(trbed: str, chrom: Optional[str] = None) -> Dict[str, TandemRepeatRegion]:
  """Load BED4+ regions and return dict mapping TR ID -> TandemRepeatRegion."""
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


def load_allosome_ploidy_bed(ploidy_bed: Optional[str] = None) -> Dict[Tuple[str, int, int], int]:
  """Load BED file with ploidy information and return dict mapping (chrom, start, end) -> ploidy."""
  if not ploidy_bed:
    logging.info('No ploidy BED file provided; assuming diploid for all regions.')
    return {}
  ploidy_dict: Dict[Tuple[str, int, int], int] = {}
  with open_bed(ploidy_bed) as fh:
    for line in fh:
      if not line.strip() or line.startswith('#'):
        continue
      cols = line.rstrip('\n').split('\t')
      if len(cols) < 5:
        logging.warning(f'Skipping malformed ploidy BED line: {line.rstrip()}')
        continue
      chrom = cols[0]
      start = int(cols[1])
      end = int(cols[2])
      ploidy = int(cols[4])
      ploidy_dict[(chrom, start, end)] = ploidy
  logging.info(
    f'Loaded ploidy information for {len(ploidy_dict)} regions from {ploidy_bed}; assuming diploid for all other regions.'
  )
  logging.debug(f'Ploidy dict contents: {ploidy_dict}')
  return ploidy_dict


def get_ploidy_for_region(ploidy_dict: Dict[Tuple[str, int, int], int], chrom: str, start: int, end: int) -> int:
  """Get ploidy for a given region from the ploidy dict, defaulting to 2 if not found."""
  if chrom not in [_[0] for _ in ploidy_dict.keys()]:
    return 2
  for (c, s, e), p in ploidy_dict.items():
    if c == chrom and start >= s and end <= e:
      return p
  return 2


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
  hap1: int,
  hap2: int,
  unphased: int,
  coverage: int,
  ploidy: int,
) -> str:
  """
  Determine dropout label given counts, chrom, sex and haplotype coverage threshold.

  return one of '', 'HaplotypeDropout', 'PhasingDropout', 'FullDropout'.
  '' indicates that both haplotypes meet coverage threshold.
  'HaplotypeDropout' indicates one haplotype below coverage threshold.
  'PhasingDropout' indicates both haplotypes below coverage threshold but total coverage sufficient.
  'FullDropout' indicates both haplotypes below coverage threshold.
  """
  if ploidy not in (0, 1, 2):
    raise ValueError(f'Unsupported ploidy value: {ploidy}')
  if ploidy == 0:
    return ''
  if ploidy == 1:
    # only one haplotype expected, check total coverage
    return 'FullDropout' if (hap1 + hap2 + unphased) < coverage * ploidy else ''
  if ploidy == 2:
    if hap1 < coverage and hap2 < coverage:
      # if both haplotypes are below coverage, call full dropout or phasing dropout depending on total coverage
      return 'FullDropout' if (hap1 + hap2 + unphased) < coverage * ploidy else 'PhasingDropout'
    if hap1 < coverage or hap2 < coverage:
      # one haplotype is not represented, call haplotype dropout
      return 'HaplotypeDropout'
    return ''
  raise ValueError(
    f'Logic error in determine_dropout: hap1={hap1}, hap2={hap2}, unphased={unphased}, coverage={coverage}, ploidy={ploidy}'
  )


def compute_and_print_results(
  tandem_repeats: Dict[str, TandemRepeatRegion],
  ploidy_dict: Dict[Tuple[str, int, int], int],
  coverage: int,
  print_all: bool,
) -> None:
  """Compute dropouts and write tab-delimited output to stdout."""
  print(
    '\t'.join(
      [
        'chrom',
        'start',
        'end',
        'trid',
        'expected_ploidy',
        'hap1_count',
        'hap2_count',
        'unphased_count',
        'fail_read_count',
        'dropout',
      ]
    ),
    file=sys.stdout,
  )
  for tr in tandem_repeats.values():
    unphased = tr.reads['0']
    hap1 = tr.reads['1']
    hap2 = tr.reads['2']
    fail = tr.reads['fail']
    ploidy = get_ploidy_for_region(ploidy_dict, tr.chrom, int(tr.chrom_start), int(tr.chrom_end))
    dropout = determine_dropout(hap1, hap2, unphased, coverage, ploidy)

    # print only regions with dropouts or fail_reads unless print_all is set
    if not print_all and dropout == '' and fail == 0:
      continue
    print(
      f'{tr.chrom}\t{tr.chrom_start}\t{tr.chrom_end}\t{tr.label}\t{ploidy}\t{hap1}\t{hap2}\t{unphased}\t{fail}\t{dropout}',
      file=sys.stdout,
    )


def main(argv=None):
  parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('trbed', help='Regions of interest, BED4+ (chrom, chromStart, chromEnd, label)')
  parser.add_argument('spanning_bam', help='Spanning reads bam.')
  parser.add_argument('--chrom', help='Chromosome to analyze', default=None)
  parser.add_argument('--ploidybed', '-p', type=str, help='BED file with ploidy information for allosomes')
  parser.add_argument('--coverage', '-c', type=int, default=2, help='Minimum per-haplotype coverage for region.')
  parser.add_argument('--print-all', '-a', action='store_true', help='Print all regions, not just those with dropouts.')
  parser.add_argument('--version', '-V', action='version', version=f'%(prog)s (version {__version__})')
  args = parser.parse_args(argv)

  logging.info(f'Starting find_trgt_dropouts (version {__version__})')

  tandem_repeats = load_trgt_catalog(args.trbed, args.chrom)
  ploidy_dict = load_allosome_ploidy_bed(args.ploidybed)
  assign_haplotags(tandem_repeats, args.spanning_bam, args.chrom)
  compute_and_print_results(tandem_repeats, ploidy_dict, args.coverage, args.print_all)
  logging.info('Completed processing')


if __name__ == '__main__':
  main()
