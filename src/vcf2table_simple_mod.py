#! /Users/karol/miniconda3/bin/python
#
#  SIMPLE MODE:
#  convert vcf to table

from multiprocessing import Pool, freeze_support

import sys, getopt, re
import pandas as pd
import numpy as np
import argparse
import vcf
import sys
import math

from datetime import datetime
from io import StringIO
from itertools import repeat

verbose = False
parent_parser = argparse.ArgumentParser()

subparsers = parent_parser.add_subparsers(title = "mode", help = "sub-command help",dest = 'mode')
subparsers.required = True

#output goes to stdout
parser_simple = subparsers.add_parser("simple", help ="run simple vcf to table conversion")
parser_simple.add_argument("-i", "--input_vcf", help = "input vcf to be parsed", required = True)
parser_simple.add_argument("-t", "--tumor", help = "ID of tumor sample", required = True)
parser_simple.add_argument("-report_normal", help = "if normal also in VCF file", action = "store_true")
parser_simple.add_argument('-debug', action = "store_true")
parser_simple.add_argument('-tab' ,action = "store_true")
parser_simple.add_argument('--build', help = "genome build")
parser_simple.add_argument('--log_time', action = "store_true")
parser_simple.add_argument('--parallel_proc', help = "Number of processes to spawn", default=1, type=int)

args = parent_parser.parse_args()

def log(msg: str):
  if(args.log_time):
    ts = datetime.now().strftime('%H:%M:%S.%f')
    print(f'{ts}: {msg}')

build = 'hg38'
build = '.'
if args.build :
    build = args.build

def eprint(str):
    print(str, file = sys.stderr)

consequence_order = \
      { "transcript_ablation":1,
        "splice_acceptor_variant":2,
        "splice_donor_variant":3,
        "stop_gained":4,
        "frameshift_variant":5,
        "stop_lost":6,
        "start_lost":7,
        "transcript_amplification":8,
        "inframe_insertion":9,
        "inframe_deletion":10,
        "missense_variant":11,
        "protein_altering_variant":12,
        "splice_region_variant":13,
        "incomplete_terminal_codon_variant":14,
        "start_retained_variant":15,
        "stop_retained_variant":16,
        "synonymous_variant":17,
        "coding_sequence_variant":18,
        "mature_miRNA_variant":19,
        "5_prime_UTR_variant":20,
        "3_prime_UTR_variant":21,
        "non_coding_transcript_exon_variant":22,
        "intron_variant":23,
        "NMD_transcript_variant":24,
        "non_coding_transcript_variant":25,
        "upstream_gene_variant":26,
        "downstream_gene_variant":27,
        "TFBS_ablation":28,
        "TFBS_amplification":29,
        "TF_binding_site_variant":30,
        "regulatory_region_ablation":31,
        "regulatory_region_amplification":32,
        "feature_elongation":33,
        "regulatory_region_variant":34,
        "feature_truncation":35,
        "intergenic_variant":36 }

def source_value(source):
    ret = 3
    if source == "RefSeq":
        ret = 2
    if source == "Ensembl":
        ret = 1
    return(ret)

def sort_annotation2(ann, found_anno, ann_header):
  ann_header.append('Variant_Classification')
  ann_header.append('ord_cs')

  rows, cols = ann.shape

  variant_classifier = lambda x: re.split('&',x)[0]
  arr = np.array([variant_classifier(xi) for xi in ann[:, ann_header.index('Consequence')]])
  ann = np.c_[ann, arr] #Variant_Classification
  con_ordering = lambda x: consequence_order[x]
  con_ord = np.array([con_ordering(xi) for xi in ann[:, ann_header.index('Variant_Classification')]], dtype=int)
  ann = np.c_[ann, con_ord] #ord_cs

  if "SOURCE" in found_anno:
    ann_header.append('ord_source')
    src_ord = np.array([source_value(xi) for xi in ann[:, ann_header.index('SOURCE')]], dtype=int)
    ann = np.c_[ann, src_ord] #ord_source
  else:
    ann = np.c_[ann, np.zeros((rows,), dtype=int)] #ord_source

  if "CANONICAL" in found_anno:
    ann_header.append('ord_canon')
    canon_ordering = lambda x: 1 if (x=="YES") else 2
    canon_ord = np.array([canon_ordering(xi) for xi in ann[:, ann_header.index('CANONICAL')]], dtype=int)
    ann = np.c_[ann, canon_ord] #ord_canon
  else:
    ann = np.c_[ann, np.zeros((rows,), dtype=int)] #ord_canon

  ann = ann[np.lexsort((ann[:, ann_header.index('ord_canon')].astype(int),
                        ann[:, ann_header.index('ord_source')].astype(int),
                        ann[:, ann_header.index('ord_cs')].astype(int)))]

  return ann, ann_header

def squish2(ann, ann_header_cp):
  log('\tSquish start')

  log('\tto string')

  with StringIO() as out:
    for item in ann[:, [ann_header_cp.index('Feature'), ann_header_cp.index('HGVSc'), ann_header_cp.index('HGVSp')]].tolist():
      out.write(';'.join(filter(None, item)) + '|')
    return out.getvalue().strip('|')

log('opening input file')
def wrap(text):
    if ',' in text:
        return "\"{}\"".format(text)
    return(text)

new_variant_I = { 'comment':'', 'gene_symbol':'.', 'genome_build':build,
                  'Chromosome':'.', 'Start_Position':'.',
                  'Variant_Classification':'.', 'Variant_Type':'.',
                  'Reference_Allele':'.', 'Tumor_Seq_Allele2':'.',
                  'HGVSc':'.', 'HGVSp':'.', 'feature_type':'.',
                  'Transcript_ID':'.', 'Exon_Number':'.'}

new_variant_II={'alt_transcripts':'.', 'Codons':'.', 'Existing_variation':'.', 'BIOTYPE':'.',
                'CANONICAL':'.',
                'SIFT':'.', 'PolyPhen':'.', 'clinvar':'.',
                'IMPACT':'.',
                'af':'.',
                'eur_af':'.', 'gnomad_af':'.', 'gnomad_nfe_af':'.',
                'max_af':'.', 'max_af_pops':'.', 'pubmed':'.'}

new_variant = {**new_variant_I,
              'tumor_DP':'.',  'tumor_AD_ref':'.',  'tumor_AD_alt':'.',  'tumor_AF':'.',
              'normal_DP':'.', 'normal_AD_ref':'.', 'normal_AD_alt':'.', 'normal_AF':'.',
              'strand_bias':'.',
              **new_variant_II}

def main():
  log('Starting')
  start = datetime.now()

  try:
    with(open(args.input_vcf, 'r')) as f:
      lines = f.readlines()
  except IOError as ioe:
    print(f'Could not read file {args.input_vcf}', file=sys.stderr)
    print(f'{ioe}', file=sys.stderr)
    sys.exit(1)

  log('file read')

  add_variant = new_variant.copy()
  print(','.join(add_variant.keys()))

  log('Starting processing')
  csq_pattern = re.compile(r"CSQ=(.+?);")
  clinvar_pattern = re.compile(r"clinvar_sig=(.+?);")
  saaf_pattern = re.compile(r"SAAF=(.+?);")
  sapp_pattern = re.compile(r"SAPP=(.+?);")

  line = lines[0].strip()
  header_lines_count = 0
  ann_header = found_anno = n_pos = t_pos = None
  while line.startswith('#'):
    ann_header, found_anno, n_pos, t_pos = process_header_line(ann_header, found_anno, line.strip())
    line = lines[header_lines_count]
    header_lines_count += 1

  lines = lines[header_lines_count - 1:]

  if args.parallel_proc == 1:
    for line in lines:
      out = process_line(line, ann_header, found_anno, t_pos, csq_pattern, clinvar_pattern, saaf_pattern, sapp_pattern)
      print(out)
  elif args.parallel_proc > 1:
    numthreads = args.parallel_proc

    pool = Pool(processes=numthreads)
    numlines = math.floor(len(lines) / numthreads)

    x = range(0, len(lines), numlines)
    params = ((ann_header, found_anno, t_pos, csq_pattern, clinvar_pattern, saaf_pattern, sapp_pattern))
    result_list = pool.starmap(worker, zip(x, (lines[line:line + numlines] for line in range(0, len(lines), numlines)), repeat(params)))

    result = {}
    for d in result_list:
      result.update(d)

    for i in result.values():
      print(i)
  else:
    print(f'Invalid number of parallel processes given: {args.parallel_proc}')
    sys.exit(1)

  end = datetime.now()
  print(f'took: {(end-start).total_seconds() * 1000} ms', file=sys.stderr)

def worker(x, lines, args):
    result = {}
    i = 0
    # print(f'worker: {x}, processing {len(lines)} lines')
    for line in lines:
      p = list(args)
      p.insert(0, line)
      result[x + i] = process_line(*p)
      i += 1
    return result

def process_line(line, ann_header, found_anno, t_pos, csq_pattern, clinvar_pattern, saaf_pattern, sapp_pattern):
  log('processing line')
  line = line.strip()
  add_variant = new_variant.copy()
  # empty new variant dictionary
  for key in add_variant:
    add_variant[key] = '.'
  add_variant['genome_build'] = build
  add_variant['comment'] = ''

  log('split 1 start')
  fields = re.split(r'\t+', line)
  log('split 1 end')

  INFO_field = fields[7] + ";"
  chr = fields[0]
  pos = fields[1]
  # eprint("chr: {}, pos: {}.".format(chr,pos))
  ref_base = fields[3]
  alt_base = fields[4]
  add_variant['Chromosome'] = chr
  add_variant['Start_Position'] = pos
  add_variant['Reference_Allele'] = ref_base
  add_variant['Tumor_Seq_Allele2'] = alt_base

  log('split 2 start')
  gth = re.split(r':', fields[8])
  log('split 2 end')

  log('simple line processing start')
  gtf = re.split(r':', fields[t_pos])

  ad_pos = gth.index("AD")
  dp_pos = gth.index("DP")

  if "AF" in gth:
    af_pos = gth.index("AF")
    af_val = gtf[af_pos]
  else:
    af_val = "."
    if gtf[ad_pos] != "." and len(re.split(',', gtf[ad_pos])) > 1 and gtf[dp_pos] != "." and gtf[dp_pos] != "0":
      af_val = "{0:.3f}".format(float(re.split(',', gtf[ad_pos])[1]) / float(gtf[dp_pos]))

  add_variant['tumor_AD_ref'] = re.split(',', gtf[ad_pos])[0]
  add_variant['tumor_AD_alt'] = re.split(',', gtf[ad_pos])[1]
  add_variant['tumor_DP'] = gtf[dp_pos]
  add_variant['tumor_AF'] = af_val
  log('simple line processing end')

  if (args.report_normal):
    gtf = re.split(r':', fields[n_pos])
    if "AF" in gth:
      af_pos = gth.index("AF")
      af_val = gtf[af_pos]
    else:
      af_val = "."
      if gtf[ad_pos] != "." and len(re.split(',', gtf[ad_pos])) > 1 and gtf[dp_pos] != "." and gtf[dp_pos] != "0":
        af_val = "{0:.3f}".format(float(re.split(',', gtf[ad_pos])[1]) / float(gtf[dp_pos]))

    add_variant['normal_AD_ref'] = re.split(',', gtf[ad_pos])[0]
    add_variant['normal_AD_alt'] = re.split(',', gtf[ad_pos])[1]
    add_variant['normal_DP'] = gtf[dp_pos]
    add_variant['normal_AF'] = af_val

  log('starting anno')

  ### Annotation
  search = re.search(csq_pattern, INFO_field)
  log('searched')
  if search != None:
    annotation = search.group(1)
    annotations = annotation.split(',')
    ann_list = []

    log('splitting')
    for a in annotations:
      spl = re.split("\|", a)
      ann_list = ann_list + [spl]

    if ann_header != "":
      log('creating numpy array')
      ann_np = np.array(ann_list)

      log('sorting')
      ann_header_cp = ann_header.copy()
      sel2, ann_header_cp = sort_annotation2(ann_np, found_anno, ann_header_cp)

    log('iterating annos')

    for anno in found_anno:
      # eprint(anno)
      if anno == "HGVSc" or anno == "HGVSp":
        coding = sel2[0][ann_header_cp.index(anno)]
        if coding != "":
          coding = re.split(":", coding)[1]
          add_variant[anno] = coding.replace('%3D', '=')
        continue

      if anno in {"MAX_AF", "gnomAD_AF", "gnomAD_NFE_AF"}:
        af = sel2[0][ann_header_cp.index(anno)]
        if af != "":
          add_variant[anno.lower()] = str("%f" % float(af))
        continue

      if anno in {"PUBMED", "MAX_AF_POPS", "AF", "EUR_AF", "Feature_type"}:
        add_variant[anno.lower()] = sel2[0][ann_header_cp.index(anno)]
        continue

      if anno == "SYMBOL":
        add_variant["gene_symbol"] = sel2[0][ann_header_cp.index(anno)]
        continue
      if anno == "Feature":
        add_variant["Transcript_ID"] = sel2[0][ann_header_cp.index(anno)]
        continue
      if anno == "EXON":
        add_variant["Exon_Number"] = sel2[0][ann_header_cp.index(anno)]
        continue
      if anno == "VARIANT_CLASS":
        add_variant['Variant_Type'] = sel2[0][ann_header_cp.index(anno)]
        continue
      if anno == "Consequence":
        add_variant["Variant_Classification"] = sel2[0][ann_header_cp.index(anno)]
        continue
      if anno == 'CLIN_SIG':
        add_variant["clinvar"] = sel2[0][ann_header_cp.index(anno)]
        continue
      if anno == 'SOURCE':
        continue
      else:
        add_variant[anno] = sel2[0][ann_header_cp.index(anno)]

      log('squish and replace')

      add_variant['alt_transcripts'] = squish2(ann_np, ann_header_cp).replace('%3D', '=')

  log('starting clinvar')

  clinvar = re.search(clinvar_pattern, INFO_field)
  if clinvar != None:
    add_variant["clinvar"] = clinvar.group(1)

  log('starting strand bias')

  # strand bias
  saaf = re.search(saaf_pattern, INFO_field)
  if saaf != None:
    sapp = re.search(sapp_pattern, INFO_field)
    if sapp != None:
      add_variant['strand_bias'] = "\"{};{}\"".format(sapp.group(1), saaf.group(1))

  add_variant.pop('SOURCE', None)
  log('line processed')
  return   ','.join(add_variant.values())

def process_header_line(ann_header, found_anno, line):
  if ("ID=CSQ" in line):
    m = re.search("Format: (.+?)\">", line)
    if m:
      ann_header_s = m.group(1)
    else:
      eprint("WARNING! Unable to extract VEP annotatin found!")
      ann_header_s = ""

    ann_header = re.split("\|", ann_header_s)
    expected_anno = ["SYMBOL", "HGVSc", "HGVSp", "IMPACT", "BIOTYPE", "CANONICAL", "Feature_type", "Feature", \
                     "EXON", "Codons", "Consequence", "SIFT", "PolyPhen", "PUBMED", "Existing_variation", 'MAX_AF_POPS', \
                     'VARIANT_CLASS', 'AF', 'EUR_AF', 'MAX_AF', 'gnomAD_AF', 'gnomAD_NFE_AF', "SOURCE", "CANONICAL",
                     "CLIN_SIG"]

    found_anno = expected_anno.copy()
    for anno in expected_anno:
      if anno not in ann_header:
        found_anno.remove(anno)
        eprint("WARNING! '{}' missing from VEP annotation!".format(anno))

    # found_anno.sort(key = lambda x: ann_header.index(x))

  n_pos = None
  t_pos = None

  if (line.startswith("#CHROM")):
    if (ann_header == ""):
      eprint("WARNING! No VEP annotatin found!")

    fields = re.split(r'\t', line)

    if (not args.tumor in fields):
      eprint("" + args.tumor + " not in sample list.")
      eprint(fields)
      exit(1)
    t_pos = fields.index(args.tumor)
    n_pos = 9 + (10 - t_pos)

  return ann_header, found_anno, n_pos, t_pos


if __name__ ==  "__main__":
  freeze_support()
  main()