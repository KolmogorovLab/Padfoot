#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 14:30:38 2025

@author: keskusa2
"""

import pysam
from collections import defaultdict
import subprocess
from Bio import Align
import os
import logging

from padfoot.annotation import load_annotation
from padfoot.copynumber import (read_wakhan_profile, gene_cn, gene_cn_fields, coverage_gaps,
                                GENE_CN_COLUMNS, GAIN_HIGH_MULTIPLIER)

logger = logging.getLogger()

class SV(object):
    __slots__ = ("bp_1", "direction_1", "bp_2", "direction_2", "supp",'supp_read_ids', 'has_ins', 'sv_type','vaf','loh', 'prec',\
                 'vcf_id', 'is_single', 'vntr', 'cluster_id', 'detailed_type', 'ins_seq', 'genes', 'repeat', 'microh', 'ins_align',\
                     'telomere', 'repeat_bp', 'hp1', 'hp2', 'cancer', 'impact',\
                     'bp_hits', 'primary', 'svlen', 'fusion', 'pairs')
    def __init__(self, bp_1, direction_1, bp_2, direction_2, supp, has_ins, sv_type,vaf, vcf_id, is_single, vntr, ins_seq, cluster_id, detailed_type, hp1, hp2, svlen=None):
        self.bp_1 = bp_1
        self.direction_1 = direction_1
        self.bp_2 = bp_2
        self.direction_2 = direction_2
        self.supp = supp
        self.supp_read_ids = ''
        self.ins_seq = ins_seq
        self.has_ins = has_ins
        self.hp1 = hp1
        self.hp2 = hp2
        self.vaf = vaf
        self.loh = ''
        self.prec = ''
        self.vcf_id = vcf_id
        self.sv_type = sv_type
        self.is_single = is_single
        self.vntr = vntr
        self.cluster_id = cluster_id
        self.detailed_type = detailed_type
        self.genes = []
        self.bp_hits = [[], []]
        self.primary = [None, None]
        self.svlen = svlen
        self.fusion = ''
        self.pairs = []
        self.repeat = []
        self.microh = []
        self.ins_align = []
        self.telomere = ''
        self.repeat_bp = ['', '']
        self.cancer = [' ',' ']
        self.impact = ''
        
    def get_jun(self):
         f1 = 'T' if self.direction_1 == '-' else 'H'
         f2 = 'T' if self.direction_2 == '-' else 'H'
         return f1+f2
     
    def get_name(self):
        label_1 = "{0}:{1}".format(self.bp_1[0], self.bp_1[1])
        label_2 = "{0}:{1}".format(self.bp_2[0], self.bp_2[1])
        return label_1 + "|" + label_2
    
    def to_str(self):
        gg = self.genes[0] if self.genes[0] else [' ', ' ', ' ']
        gene1 = '\t'.join(gg)
        
        gg = self.genes[1] if self.genes[1] else [' ', ' ', ' ']
        gene2 = '\t'.join(gg)
        
        typ = self.genes[2] if len(self.genes)==3 else ' '
        microh = self.microh[0] if self.microh else ' '
        tel = str(self.telomere) if self.telomere else ' '
        ins_al = ';'.join(set(self.ins_align)) if self.ins_align else ' '
        rep = ','.join([':'.join([s[0], str(s[1])]) for s in self.repeat])

        st = self.direction_1 + self.direction_2
        return '\t'.join([self.vcf_id, self.bp_1[0], str(self.bp_1[1]), self.bp_2[0], str(self.bp_2[1]), str(self.supp),\
         str(round(self.vaf, 3)), st, self.impact or ' ', self.cancer[0], gene1,  self.repeat_bp[0], self.cancer[1], gene2, self.repeat_bp[1], typ, microh, str(self.vntr), tel, ins_al, rep]
         + self.new_fields())

    def other_genes(self, i):
        """Genes overlapping breakend i (body or padding) other than its primary gene."""
        primary = self.primary[i].gene.key if self.primary[i] else None
        return ','.join(dict.fromkeys(h.gene.name for h in self.bp_hits[i] if h.gene.key != primary))

    def new_fields(self):
        whole = [p.gene.name for p in self.pairs if any(l in WHOLE_GENE_LABELS.values() for l in p.labels)]
        whole = WHOLE_GENE_LABELS[self.sv_type] + ':' + ','.join(whole) if whole else ''
        alt = ','.join(dict.fromkeys(p.gene.name for p in self.pairs if 'alt_isoform_exon' in p.labels))
        known = 'yes' if len(self.genes) == 3 and self.genes[2] == 'oncogenic_fusion' else 'no'
        return [self.other_genes(0), self.other_genes(1), self.fusion, known, whole, alt,
                str(self.cluster_id or ''), str(self.detailed_type or ''), self.haplotype()]

    def haplotype(self):
        """HP1/HP2 from the phased breakend(s); HP0 if neither is phased; HP1|HP2 if they disagree."""
        hps = sorted({self.hp1, self.hp2} - {0})
        return '|'.join('HP%d' % h for h in hps) if hps else 'HP0'

    def pos_str(self):
        if self.bp_1[0] == self.bp_2[0]:
            return self.bp_1[0] + ':' + str(self.bp_1[1]) + '-' + str(self.bp_2[1])
        return self.bp_1[0] + ':' + str(self.bp_1[1]) + '-' + self.bp_2[0] + ':' + str(self.bp_2[1])

class Gene(object):
    __slots__ = ("ref_id", "pos_1", "pos_2", "Gene_symbol", "gene_id", "cn", 'pairs', 'cancer')
    def __init__(self, ref_id, pos_1, pos_2, Gene_symbol, gene_id=''):
        self.gene_id = gene_id
        self.ref_id = ref_id
        self.pos_1 = pos_1
        self.pos_2 = pos_2
        self.Gene_symbol = Gene_symbol
        self.cn = None
        self.pairs = []
        self.cancer = ''
    def impact(self, hp):
        """Most severe SV tier on one haplotype (0 = unphased)."""
        return most_severe(p.tier for p in self.pairs if hp in p.haplotypes)
    def to_str(self):
        pos = self.ref_id + ':' + str(self.pos_1) + '-' + str(self.pos_2)
        sv_info = []
        for hp in (0, 1, 2):
            svs = sorted(set(p.to_sum() for p in self.pairs if hp in p.haplotypes))
            sv_info.append(','.join(svs))
        svs = {p.sv.vcf_id: p.sv for p in self.pairs}
        in_cluster = 'yes' if any(sv.cluster_id for sv in svs.values()) else 'no'
        return '\t'.join([self.Gene_symbol, pos, self.cancer, self.impact(1), sv_info[1], self.impact(2), sv_info[2],
                          self.impact(0), sv_info[0], self.gene_id, str(len(svs)), in_cluster] + gene_cn_fields(self.cn))


class GeneSVPair(object):
    """One SV as seen from one gene: the class and tier for this gene."""
    __slots__ = ("gene", "sv", "hits", "cls", "tier", "labels")
    def __init__(self, gene, sv, hits, cls, tier, labels=()):
        self.gene = gene
        self.sv = sv
        self.hits = hits                # (hit at bp1 or None, hit at bp2 or None)
        self.cls = cls
        self.tier = tier
        self.labels = list(labels)
    @property
    def haplotypes(self):
        hps = {hp for hit, hp in zip(self.hits, (self.sv.hp1, self.sv.hp2)) if hit}
        return hps or {self.sv.hp1, self.sv.hp2}
    def to_sum(self):
        fus = ',' + self.sv.fusion if self.cls in FUSION_CLASSES else ''
        what = self.cls or ','.join(self.labels)
        return (self.tier or '-') + ':[' + self.sv.vcf_id + ',' + self.sv.pos_str() + fus + ',' + what + ']'

def get_bps(vcf_file):
    svs = defaultdict(list)
    vcf = pysam.VariantFile(vcf_file)
    for var in vcf:
        if var.info['SVTYPE'] == 'BND' and var.id.endswith('_2'):
            continue    # the _1 record carries both breakends (HP=<bp1>|<bp2>)
        bp_1 = (var.chrom, var.pos)
        dir_ls= var.info['STRANDS'] if 'STRANDS' in var.info.keys() else ('','')
        dir1, dir2 = (dir_ls[0], dir_ls[-1])
        
        ins_seq = '' if not var.info['SVTYPE'] == 'INS'  else var.alts[0]
        sample_id = var.samples.keys()[0]
        supp = var.samples[sample_id]['DV']
        vaf = var.samples[sample_id]['VAF']
        vntr = True if 'INSIDE_VNTR' in var.info.keys() else False
        has_ins = '' if not 'INSSEQ' in var.info.keys() else var.info['INSSEQ']
        is_single = False if not var.info['SVTYPE'] == 'sBND' else True
        sv_type = var.info['SVTYPE']
        vcf_id = var.id
        hp1 = hp2 = 0
        bp_2 = bp_1
        if 'HP' in var.info.keys():
            hp1, hp2 = [int(h) for h in var.info['HP'].split('|')]
        
        if var.info['SVTYPE'] == 'sBND':
            chr2 = var.chrom
            
        elif var.info['SVTYPE'] == 'BND':
            dir1, dir2 = var.info['STRANDS']
            chr2, pos2 = var.alts[0].replace('[','').replace(']','').replace('N', '').split(':')
            pos2 = int(pos2)
            bp_2 = (chr2, pos2)
            
        elif var.info['SVTYPE'] in ['DEL', 'DUP', 'INV']:
            end_pos = var.stop if var.stop else var.pos + var.info['SVLEN']
            bp_2 = (var.chrom, end_pos)
            
        elif var.info['SVTYPE'] == 'INS':
            bp_2 = bp_1
            
        cluster_id = ''
        if 'CLUSTERID' in var.info.keys():
            cluster_id = var.info['CLUSTERID']
            
        detailed_type = ''
        if 'DETAILED_TYPE' in var.info.keys():
            detailed_type = var.info['DETAILED_TYPE']
            
        svlen = var.info['SVLEN'] if 'SVLEN' in var.info.keys() else None
        svs[vcf_id]=(SV(bp_1, dir1, bp_2, dir2, supp, has_ins, sv_type,vaf, vcf_id, is_single, vntr,ins_seq, cluster_id, detailed_type, hp1, hp2, svlen))
    return list(svs.values())

def get_genes(annotation_file):
    """Load the compact gene annotation (see padfoot/annotation.py)."""
    return load_annotation(annotation_file)

def pick_primary(hits, other_hits):
    """Primary gene for one breakend, in order of preference:
    1. a gene whose body contains the breakend (over one reached only through the padding);
    2. a gene whose body also contains the other breakend, then one the other
       breakend reaches at all (so that both ends agree where possible);
    3. the shortest gene (a nested gene over its host)."""
    if not hits:
        return None
    other_body = {h.gene.key for h in other_hits if h.in_body}
    other_any = {h.gene.key for h in other_hits}
    return min(hits, key=lambda h: (not h.in_body, h.gene.key not in other_body, h.gene.key not in other_any,
                                    h.gene.end - h.gene.start, h.gene.start, h.gene.gene_id))

TIER_RANK = {'HIGH': 2, 'MODERATE': 1, '': 0}
FUSION_CLASSES = ('possible_fusion', 'oncogenic_fusion')
WHOLE_GENE_LABELS = {'DEL': 'whole_gene_deletion', 'DUP': 'whole_gene_duplication', 'INV': 'whole_gene_inversion'}
WHOLE_GENE_MAX_SV = 1_000_000

def most_severe(tiers):
    return max(tiers, key=TIER_RANK.get, default='')

def is_5prime(hit, direction):
    """True if the breakend keeps the 5' part of the gene. Severus strand '+'
    keeps the sequence to the left of the breakend, '-' the sequence to its right."""
    return (direction == '+') == (hit.gene.strand == '+')

def sv_length(sv):
    if sv.sv_type == 'INS':
        return sv.svlen or 0
    if sv.bp_1[0] != sv.bp_2[0]:
        return 0
    return abs(sv.bp_2[1] - sv.bp_1[1])

def within_gene_class(sv, h1, h2):
    """Class of an SV with both breakends in (or padding of) the same gene."""
    if h1.key == h2.key:
        if h1.region == 'exon':
            return 'frameshift' if sv_length(sv) % 3 else 'exon_altering'
        if h1.region == 'intron':
            return 'intronic'
        return 'promoter/UTR'
    return 'between_exons'

def sv_class(sv, p1, p2):
    """SV-level class from the primary gene at each breakend (v0.1 classes)."""
    if not p1 and not p2:
        return ''
    if not p1 or not p2:
        return 'Gene-NonCoding'
    if p1.gene.key == p2.gene.key:
        return within_gene_class(sv, p1, p2)
    if p1.in_body and p2.in_body and is_5prime(p1, sv.direction_1) != is_5prime(p2, sv.direction_2):
        return 'possible_fusion'
    return 'between_genes'

def pair_tier(sv, cls, h1, h2):
    """Tier of an SV for one gene (class x location x SV type); never cancer-dependent."""
    if cls in ('exon_altering', 'frameshift'):
        return 'HIGH' if h1.coding else 'MODERATE'
    if cls == 'between_exons':
        if h1.region == 'flank' and h2.region == 'flank':
            return ''           # the whole gene lies between the breakends
        return 'HIGH'
    if cls in FUSION_CLASSES:
        return 'HIGH'
    if cls in ('intronic', 'promoter/UTR'):
        return 'MODERATE'
    hit = h1 or h2              # between_genes / Gene-NonCoding: one end in this gene
    if not hit.in_body or sv.sv_type == 'DUP':
        return 'MODERATE'
    return 'HIGH'

def gene_sv_pairs(annotation, sv, p1, p2, whole_gene_max):
    hits = [{h.gene.key: h for h in sv.bp_hits[0]}, {h.gene.key: h for h in sv.bp_hits[1]}]
    primary = {p.gene.key for p in (p1, p2) if p}
    genes = {}
    for hs in hits:
        for key, h in hs.items():
            genes[key] = h.gene
    whole = set()
    if (sv.sv_type in WHOLE_GENE_LABELS and sv.bp_1[0] == sv.bp_2[0]
            and sv_length(sv) <= whole_gene_max):
        for g in annotation.genes_within(sv.bp_1[0], min(sv.bp_1[1], sv.bp_2[1]), max(sv.bp_1[1], sv.bp_2[1])):
            whole.add(g.key)
            genes.setdefault(g.key, g)

    pairs = []
    for key, gene in sorted(genes.items(), key=lambda kv: (kv[1].start, kv[1].end, kv[1].gene_id)):
        h1, h2 = hits[0].get(key), hits[1].get(key)
        if h1 and h2:
            cls = within_gene_class(sv, h1, h2)
        elif h1 or h2:
            other = p2 if h1 else p1
            if other is None:
                cls = 'Gene-NonCoding'
            elif sv.genes[2] in FUSION_CLASSES and key in primary:
                cls = sv.genes[2]
            else:
                cls = 'between_genes'
        else:
            cls = ''
        labels = []
        if key in whole:
            labels.append(WHOLE_GENE_LABELS[sv.sv_type])
            tier = ''
        else:
            tier = pair_tier(sv, cls, h1, h2)
        if any(h and h.alt_isoform_exon for h in (h1, h2)):
            labels.append('alt_isoform_exon')
        pairs.append(GeneSVPair(gene, sv, (h1, h2), cls, tier, labels))
    return pairs

def annot_SVS(annotation, svs, by_gene, whole_gene_max=WHOLE_GENE_MAX_SV):
    for sv in svs:
        hits1 = annotation.lookup(*sv.bp_1)
        hits2 = annotation.lookup(*sv.bp_2)
        p1, p2 = pick_primary(hits1, hits2), pick_primary(hits2, hits1)
        sv.bp_hits = [hits1, hits2]
        sv.primary = [p1, p2]
        sv.genes = [p1.as_tuple() if p1 else (), p2.as_tuple() if p2 else (), sv_class(sv, p1, p2)]
        if sv.genes[2] in FUSION_CLASSES:
            five, three = (p1, p2) if is_5prime(p1, sv.direction_1) else (p2, p1)
            sv.fusion = five.gene.name + '::' + three.gene.name
        sv.pairs = gene_sv_pairs(annotation, sv, p1, p2, whole_gene_max)
        sv.impact = most_severe(p.tier for p in sv.pairs)
        for pair in sv.pairs:
            add_gene(by_gene, pair.gene).pairs.append(pair)

def add_gene(by_gene, gene):
    """by_gene row of an annotation gene, keyed by (chrom, gene_id)."""
    if gene.key not in by_gene:
        by_gene[gene.key] = Gene(gene.chrom, gene.start, gene.end, gene.name, gene.gene_id)
    return by_gene[gene.key]

def annot_gene_cn(annotation, profile, by_gene, sex='unknown', gain_high_multiplier=GAIN_HIGH_MULTIPLIER):
    """CN fields for every annotated gene (all genes get a by_gene row)."""
    for g in annotation.genes:
        add_gene(by_gene, g).cn = gene_cn(profile, g.chrom, g.start, g.end, sex, gain_high_multiplier)

def check_cn_coverage(profile, annotation, ref):
    """Warn about chromosomes and regions without a CN call."""
    lengths = {}
    if os.path.exists(ref + '.fai'):
        with open(ref + '.fai') as f:
            for line in f:
                name, length = line.split('\t')[:2]
                lengths[name] = int(length)
    chroms = list(dict.fromkeys(list(profile.by_chrom) + [g.chrom for g in annotation.genes]))
    absent, gaps = coverage_gaps(profile, lengths, chroms)
    if absent:
        logger.warning("No CN profile for %d chromosome(s) with annotated genes: %s", len(absent), ', '.join(absent))
    if gaps:
        total = sum(e - s + 1 for _, s, e in gaps)
        logger.warning("CN profile has %d gap(s) totalling %d bp (no call for both haplotypes); first: %s",
                       len(gaps), total, ', '.join(f'{c}:{s}-{e}' for c, s, e in gaps[:5]))
    return absent, gaps

def run_command(cmd):
    p = subprocess.Popen(cmd, shell= True)
    p.wait()

def write_ins(svs, ref, t, repeatmasker_species=None):
    fa_out = open('temp_ins.fa', 'w')
    for sv in svs:
        if sv.ins_seq:
            fa_out.write('>' + sv.vcf_id + '\n')
            fa_out.write(sv.ins_seq)
            fa_out.write('\n')
        elif sv.has_ins:
            fa_out.write('>' + sv.vcf_id + '\n')
            fa_out.write(sv.has_ins)
            fa_out.write('\n')
    fa_out.close()
    if os.path.exists('temp_ins.fa.out'):
        os.remove('temp_ins.fa.out')            # stale output of an earlier run in the same temp dir
    if repeatmasker_species:
        run_command(f"RepeatMasker -species {repeatmasker_species} temp_ins.fa")
    run_command(f"minimap2 -ax map-ont {ref} temp_ins.fa -k 17 -y -K 5G -t {t} --eqx | samtools sort -@ {t} -m 4G > temp_ins.bam")
    run_command(f"samtools index -@ {t} temp_ins.bam")
    

def get_repeat(svls):
    if os.path.exists('temp_ins.fa.out'):
        rep_file = open('temp_ins.fa.out', 'r')
    else:
        logger.warning("RepeatMasker output not found (temp_ins.fa.out); RepeatMasker may have failed")
        return ''
    reps = defaultdict(int)
    svrep = defaultdict(list)
    sv_id = ''
    ll = 0
    for line in rep_file:
        ll+=1
        if ll < 4:
            continue
        l = line.split()
        if not sv_id:
            sv_id = l[4]
        if sv_id == l[4]:
            reps[l[10]]+= int(l[6]) - int(l[5])
            sv_len = int(l[7][1:-1]) + int(l[6])
        else:
            max_rep = max(reps, key=reps.get)
            replen = reps[max_rep]
            if replen > sv_len * 0.8:
                svrep[sv_id] = (max_rep, replen)
            reps = defaultdict(int)
            sv_id = l[4]
            reps[l[10]]+= int(l[6]) - int(l[5])
            sv_len = int(l[7][1:-1]) + int(l[6])
    for key,val in svrep.items():
        svls[key].repeat.append(val)

def check_homology(seq, sv):
    seq1, seq2 = seq
    THR = 25
    MINSCORE = 3
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.gap_score = -100
    aligner.mismatch_score = -100
    ss = []
    for i in [THR-10,THR]:
        for j in [THR-10,THR]:
            a = aligner.align(seq1[i:i+11], seq2[j:j+11])
            if not a:
                continue
            else:
                aa = a[0]
                if aa.score < MINSCORE:
                    continue
                pos1 = abs(THR-i-aa.aligned[0][0][0]) + abs(THR-j - aa.aligned[1][0][0])
                if aa.score-pos1 > 0: 
                    ss.append((aa.score-pos1, aa.sequences[0][aa.aligned[0][0][0]:aa.aligned[0][0][1]]))
    if ss:
        ss.sort(key=lambda s:s[0])
        sv.microh.append(ss[-1][1])
    
def get_microhomology(svs, ref):
    """Microhomology at the junction, from +-25 bp of reference around each breakend."""
    THR = 25
    missing = set()
    with pysam.FastaFile(ref) as fasta:
        for sv in svs:
            if sv.sv_type == 'INS':
                continue
            seq = []
            for chrom, pos in (sv.bp_1, sv.bp_2):
                if chrom not in fasta:
                    missing.add(chrom)
                    break
                seq.append(fasta.fetch(chrom, max(0, pos - THR - 1), pos + THR).upper())
            else:
                check_homology(seq, sv)
    if missing:
        logger.warning("Microhomology not computed for SVs on %d contig(s) missing from the reference: %s",
                       len(missing), ', '.join(sorted(missing)))

def get_align(svls):
    MAPQ_THR = 45
    aln_file = pysam.AlignmentFile('temp_ins.bam', "rb")
    for aln in aln_file.fetch():
        if not aln.is_secondary and not aln.is_unmapped and aln.mapq > MAPQ_THR:
            pos = aln.reference_name +':' + str(aln.reference_start) +'-' + str(aln.reference_end)
            svls[aln.query_name].ins_align.append(pos)
    
def get_tel(svs):
    telseq = ['CCCTAA', 'TTAGGG']
    MIN_THR = 4
    for sv in svs:
        if sv.ins_seq or sv.has_ins:
            seq = sv.ins_seq if sv.ins_seq else sv.has_ins
            c1 = seq.count(telseq[0])
            if c1 < MIN_THR:
                c1 = seq.count(telseq[1])
            if c1 >= MIN_THR:
                sv.telomere = c1

def annot_bp_repeat(svls, rm_bed):
    f = open('temp_bps.bed', 'w')
    for sv in svls.values():
        f.write('\t' .join([sv.bp_1[0], str(max(sv.bp_1[1]-5,0)), str(sv.bp_1[1]+5), sv.vcf_id, 'BP1']))
        f.write('\n')
        if not sv.sv_type == 'INS':
            f.write('\t' .join([sv.bp_2[0], str(max(0, sv.bp_2[1]-5)), str(sv.bp_2[1]+5), sv.vcf_id, 'BP2']))
            f.write('\n')
            
    run_command(f'bedtools intersect -a temp_bps.bed -b {rm_bed} -wb > temp_int.bed')
    
    f1 = open('temp_int.bed')
    for line in f1:
        l = line.split()
        ind = 0 if l[4] == 'BP1' else 1
        svls[l[3]].repeat_bp[ind] = l[-1]
        
def annot_ins(svs, ref,t, rm_bed, repeatmasker_species=None):
    svls = defaultdict(list)
    for sv in svs:
         svls[sv.vcf_id] = sv
    write_ins(svs, ref, t, repeatmasker_species)
    if repeatmasker_species:
        get_repeat(svls)
    get_align(svls)
    get_tel(svs)
    annot_bp_repeat(svls, rm_bed)
    
    
CANCER_GENE_COLUMNS = ('Gene_symbol', 'Role_in_cancer', 'fusion')

def load_cancer_genes(path):
    """Cancer gene list: a TSV with a header containing Gene_symbol, Role_in_cancer
    and fusion (comma-separated partner symbols); other columns are ignored.
    Returns ({gene: role}, {gene: [fusion partners]})."""
    roles, fusion = {}, defaultdict(list)
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        missing = [c for c in CANCER_GENE_COLUMNS if c not in header]
        if missing:
            raise ValueError(f"{path}: cancer gene list lacks column(s) {missing}")
        col = {c: header.index(c) for c in CANCER_GENE_COLUMNS}
        for line in f:
            row = line.rstrip('\n').split('\t')
            row += [''] * (len(header) - len(row))
            gene = row[col['Gene_symbol']].strip()
            if not gene:
                continue
            role = row[col['Role_in_cancer']].strip()
            if role:
                roles[gene] = role
            partners = [x.strip() for x in row[col['fusion']].split(',') if x.strip()]
            fusion[gene].extend(x for x in partners if x not in fusion[gene])
    return roles, fusion

def is_known_fusion(a, b, fusion):
    """Gene pair listed as a fusion, under either partner."""
    return b in fusion.get(a, ()) or a in fusion.get(b, ())

def cancer_annot_svs(svs, cancer_genes, fusion):
    """Cancer-gene role of each breakend's gene and the known-fusion lookup (labels only):
    a possible_fusion of a listed pair becomes oncogenic_fusion (same HIGH tier)."""
    for sv in svs:
        if not sv.genes[2]:
            continue
        for i in (0, 1):
            if sv.genes[i] and sv.genes[i][0] in cancer_genes:
                sv.cancer[i] = cancer_genes[sv.genes[i][0]]
        if sv.genes[2] == 'possible_fusion' and is_known_fusion(sv.genes[0][0], sv.genes[1][0], fusion):
            sv.genes[2] = 'oncogenic_fusion'
            for pair in sv.pairs:
                if pair.cls == 'possible_fusion':
                    pair.cls = 'oncogenic_fusion'

def cancer_annot_genes(by_gene, cancer_genes):
    for gs in by_gene.values():
        if gs.Gene_symbol in cancer_genes:
            gs.cancer = cancer_genes[gs.Gene_symbol]

def cancer_annot(svs, by_gene, cancer_gene_file):
    cancer_genes, fusion = load_cancer_genes(cancer_gene_file)
    logger.info("Cancer gene list: %s (%d genes with a role, %d with fusion partners)",
                cancer_gene_file, len(cancer_genes), sum(1 for v in fusion.values() if v))
    cancer_annot_svs(svs, cancer_genes, fusion)
    cancer_annot_genes(by_gene, cancer_genes)

SV_COLUMNS = ['SV ID', 'chr_1', 'pos_1', 'chr_2', 'pos_2', 'N_supp', 'Vaf', 'Strand', 'Impact',
              'Role in Cancer1', 'Gene Name1', 'Pos_1', 'Strand_1', 'Repeat_Anno1',
              'Role in Cancer2', 'Gene Name2', 'Pos_2', 'Strand_2', 'Repeat_Anno2',
              'Type', 'Microhomology', 'VNTR', 'Telomere_repeat', 'Aligned_pos(INS)', 'Repeat_annot',
              'Other_genes1', 'Other_genes2', 'fusion_partners', 'known_fusion', 'whole_gene_labels',
              'alt_isoform_exon', 'cluster_id', 'detailed_type', 'haplotype']

GENE_COLUMNS = ['Gene Name', 'Gene Pos', 'Role in Cancer', 'HP1 Impact', 'HP1 SV', 'HP2 Impact', 'HP2 SV',
                'Unphased annot', 'Unphased SV', 'gene_id', 'n_SVs', 'in_cluster'] + GENE_CN_COLUMNS

def output_svs(svs, out_dir):
    with open(os.path.join(out_dir, 'annotated_svs.tsv'), 'w') as f:
        f.write('\t'.join(SV_COLUMNS) + '\n')
        for sv in svs:
            f.write(sv.to_str() + '\n')

def output_genes(by_gene, out_dir):
    with open(os.path.join(out_dir, 'by_gene.tsv'), 'w') as f:
        f.write('\t'.join(GENE_COLUMNS) + '\n')
        for row in by_gene.values():
            f.write(row.to_str() + '\n')


def _location(hit):
    if hit is None:
        return ''
    if hit.region == 'exon':
        return ('coding_' if hit.coding else 'utr_') + hit.label
    if hit.region == 'intron':
        return hit.label
    return hit.side + '_padding'

GENE_SV_PAIR_COLUMNS = ['gene', 'gene_id', 'gene_pos', 'sv_id', 'sv_type', 'haplotype', 'class', 'tier',
                        'bp1_location', 'bp2_location', 'fusion', 'labels']

def output_gene_sv_pairs(svs, out_dir):
    with open(os.path.join(out_dir, 'gene_sv_pairs.tsv'), 'w') as f:
        f.write('\t'.join(GENE_SV_PAIR_COLUMNS) + '\n')
        for sv in svs:
            for p in sv.pairs:
                g = p.gene
                f.write('\t'.join([g.name, g.gene_id, f'{g.chrom}:{g.start}-{g.end}', sv.vcf_id, sv.sv_type,
                                   ','.join(map(str, sorted(p.haplotypes))), p.cls, p.tier,
                                   _location(p.hits[0]), _location(p.hits[1]),
                                   sv.fusion if p.cls in FUSION_CLASSES else '', ','.join(p.labels)]) + '\n')

def annotate_things(args):
    vcf_file, t, ref, out_dir = args.vcf_file, args.threads, args.ref, args.out_dir
    by_gene = defaultdict(list)
    annotation = get_genes(args.gff_file)
    svs = get_bps(vcf_file)
    profile = None
    if getattr(args, 'wakhan_profile', None):
        profile = read_wakhan_profile(args.wakhan_profile)
        logger.info("CN profile: %s (%d fragments, modal CN %s)", profile.path, len(profile), profile.modal_cn)
        check_cn_coverage(profile, annotation, ref)
        annot_gene_cn(annotation, profile, by_gene, getattr(args, 'sex', 'unknown'),
                      getattr(args, 'gain_high_multiplier', GAIN_HIGH_MULTIPLIER))
    annot_SVS(annotation, svs, by_gene, getattr(args, 'whole_gene_max_sv', WHOLE_GENE_MAX_SV))
    if getattr(args, 'cancer_genes', None):
        cancer_annot(svs, by_gene, args.cancer_genes)
    annot_ins(svs, ref,t, args.rm_file, getattr(args, 'repeatmasker_species', None))
    get_microhomology(svs, ref)
    output_svs(svs, out_dir)
    output_genes(by_gene, out_dir)
    output_gene_sv_pairs(svs, out_dir)
    return (annotation, profile, svs, by_gene)
