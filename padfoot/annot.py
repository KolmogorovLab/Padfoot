#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 14:30:38 2025

@author: keskusa2
"""

import pysam
from collections import defaultdict, Counter
import bisect
import subprocess
import gzip
from Bio import Align
import copy
import shutil
import numpy as np

class SV(object):
    __slots__ = ("bp_1", "direction_1", "bp_2", "direction_2", "supp",'supp_read_ids', 'has_ins', 'sv_type','vaf','loh', 'prec',\
                 'vcf_id', 'is_single', 'vntr', 'cluster_id', 'detailed_type', 'ins_seq', 'genes', 'repeat', 'microh', 'ins_align',\
                     'telomere', 'repeat_bp', 'cn_assigned', 'cn_altering', 'hp1', 'hp2')
    def __init__(self, bp_1, direction_1, bp_2, direction_2, supp, has_ins, sv_type,vaf, vcf_id, is_single, vntr, ins_seq, cluster_id, detailed_type, hp1, hp2):
        self.bp_1 = bp_1
        self.bp_2 = bp_2
        self.direction_1 = direction_1
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
        self.repeat = []
        self.microh = []
        self.ins_align = []
        self.telomere = ''
        self.repeat_bp = ['', '']
        self.cn_assigned  = ''
        self.cn_altering = False
        
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
         str(round(self.vaf, 3)), st, gene1,self.repeat_bp[0], gene2, self.repeat_bp[1], typ, microh, str(self.vntr), tel, ins_al, rep])
            
class CNA(object):
    __slots__ = ("ref_id", "pos_1", "pos_2", "cn", "haplotype",'LOH', 'sv1', 'sv2', 'genes', 'dir1')
    def __init__(self, ref_id, pos_1, pos_2, cn, haplotype, LOH, sv1, sv2):
        self.ref_id = ref_id
        self.pos_1 = pos_1
        self.pos_2 = pos_2
        self.cn = cn
        self.haplotype = haplotype
        self.LOH = LOH
        self.sv1 = sv1
        self.sv2 = sv2
        self.genes = []
        self.dir1 = ''
    def get_name(self):
        return 'HP' + str(self.haplotype) +':' + self.ref_id + ':' + str(self.pos_1) + '-' + str(self.pos_2)
    
def get_bps(vcf_file):
    svs = defaultdict(list)
    vcf = pysam.VariantFile(vcf_file)
    for var in vcf:
        if "_2" in var.id:
            continue
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
        if 'HP' in var.info.keys():
            hp1 = hp2 = var.info['HP']
        
        if var.info['SVTYPE'] == 'sBND':
            chr2 = var.chrom
            
        elif var.info['SVTYPE'] == 'BND':
            if '_2' in var.id:
                if 'HP' in var.info.keys():
                    hp2 = var.info['HP']
                    svs[var.id.replace('_2', '_1')].hp2 = hp2
                continue
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
            
        svs[vcf_id]=(SV(bp_1, dir1, bp_2, dir2, supp, has_ins, sv_type,vaf, vcf_id, is_single, vntr,ins_seq, cluster_id, detailed_type, hp1, hp2))
    return list(svs.values())

def get_CNA(cna_vcf, svs):
    vcf = pysam.VariantFile(cna_vcf)
    hp1ls = defaultdict(list)
    hp2ls = defaultdict(list)
    LOH  = defaultdict(list)
    CNAs = defaultdict(list)
    ploidy = []
    cov1 = []
    
    for var in vcf:
        ref_id, pos_1, pos_2 = var.chrom, var.pos, var.stop
        hp1, hp2 = var.samples['Sample']['CN1'], var.samples['Sample']['CN2']
        if hp1:
            cov1.append( var.samples['Sample']['COV1']/hp1)
        if ref_id in hp1ls:
            if not hp1ls[ref_id][0][-1] == hp1:
                hp1ls[ref_id][0].append(hp1)
                hp1ls[ref_id][1].append(pos_1)
            if not hp2ls[ref_id][0][-1] == hp2:
                hp2ls[ref_id][0].append(hp2)
                hp2ls[ref_id][1].append(pos_1)
        else:
            hp1ls[ref_id]= [[hp1],[pos_1]]
            hp2ls[ref_id] = [[hp2],[pos_1]]
        if hp1 == 0 or hp2 == 0:
            if ref_id in LOH.keys():
                LOH[ref_id][0].append(pos_1)
                LOH[ref_id][1].append(pos_1)
            else:
                LOH[ref_id]= [[pos_1],[pos_2]]

    cn1_cov = int(np.median(cov1)) *0.75
    
    for hp, hpls in enumerate([hp1ls, hp2ls]):
        ploidy_hp = 0
        len_cn = 0
        for ref, values in hpls.items():
            cnls = values[0]
            posls = values[1]
            loh = False
            for i, (pos_1, pos_2) in enumerate(zip(posls[:-1], posls[1:])):
                cn = cnls[i]
                ls1 = [(ref_id, pos_1),(ref_id, pos_1+1),(ref_id, pos_1-1)]
                sv1 = [sv for sv in svs if sv.bp_1 in ls1 and sv.direction_1 == '-' or sv.bp_2 in ls1 and sv.direction_2 == '-']
                sv1 = sv1[0] if sv1 else ''
                ls2 = [(ref_id, pos_2),(ref_id, pos_2+1),(ref_id, pos_2-1)]
                sv2 = [sv for sv in svs if sv.bp_1 in ls2 and sv.direction_1 == '+' or sv.bp_2 in ls2 and sv.direction_2 == '+' ]
                sv2 = sv2[0] if sv2 else ''
                if ref_id in LOH.keys() and  pos_1 in LOH[ref_id][0]:
                    loh = True
                CNAs[(ref_id, hp+1)].append(CNA(ref_id, pos_1, pos_2, cn, hp+1, loh, sv1, sv2))
                len_cn += (pos_2 - pos_1)
                ploidy_hp += cn * (pos_2 - pos_1)
        ploidy.append(round(ploidy_hp/len_cn))     

    check_hp_svs(CNAs)
    for (ref,hp), cnas in CNAs.items():
        pl = ploidy[hp-1]
        for cna in cnas:
            if cna.cn > pl:
                cna.dir1 = 'AMP'
            elif cna.cn < pl:
                cna.dir1 = 'DEL'
    check_cn_altering_svs(svs, cn1_cov)
    return(CNAs, ploidy)

def check_cn_altering_svs(svs, cn1_cov):
    for sv in svs:
        sv.cn_altering = round(sv.supp/cn1_cov)

def check_hp_svs(cnas):
    multsvs = defaultdict(list)
    for cns in cnas.values():
        for cn in cns:
            if cn.sv1:
                multsvs[cn.sv1].append(cn)
            if cn.sv2:
                multsvs[cn.sv2].append(cn)
    for sv, cns in multsvs.items():
        xx = Counter([(cn.ref_id, cn.pos_1) for cn in cns])
        xx = [x for x in xx.values() if x > 1]
        
        yy = Counter([(cn.ref_id, cn.pos_2) for cn in cns])
        yy = [x for x in yy.values() if x > 1]
        
        if not xx and not yy:
            continue
        hp2 = [cn for cn in cns if cn.haplotype == 2]
        new_sv = copy.copy(sv)
        for cn in hp2:
            if cn.sv1 == sv:
                cn.sv1 = new_sv
            else:
                cn.sv2 = new_sv
    
def get_genes(gff_file):
    THR = 10000
    genes = defaultdict(list)
    exon_pos = defaultdict(list)
    fopen =gzip.open(gff_file, 'rt')
    for line in fopen:
        ref_id, gene_name, strand, typ, start, end = line.split()
        start = int(start)
        end = int(end)
        THR1, THR2 = (THR, 0) if strand == '+' else (0, THR)
        if typ == 'gene':
            if ref_id in genes.keys():
                genes[ref_id][0].append(gene_name)
                genes[ref_id][1].append(start-THR)
                genes[ref_id][2].append(end+THR)
            else:
                genes[ref_id] = [[],[],[]]
                genes[ref_id][0].append(gene_name)
                genes[ref_id][1].append(start-THR)
                genes[ref_id][2].append(end+THR)
        elif typ.startswith(gene_name) or 'codon' in typ or 'UTR' in typ:
            continue
        else:
            if gene_name in exon_pos.keys():
                exon_pos[gene_name][0].append(start)
                exon_pos[gene_name][1].append(end)
                exon_pos[gene_name][2].append(typ)
            else:
                exon_pos[gene_name] = [[],[],[],[strand]]
                exon_pos[gene_name][0].append(start)
                exon_pos[gene_name][1].append(end)
                exon_pos[gene_name][2].append(typ)
    for gene_name,exons in exon_pos.items():
        exonls = [ind for ind, typ in enumerate(exons[2]) if 'exon' in typ]
        utrs5 = [ind for ind, typ in enumerate(exons[2]) if typ == 'five_prime_UTR']
        utrs3 = [ind for ind, typ in enumerate(exons[2]) if typ == 'three_prime_UTR']
        if utrs5:
            st5 = min([exons[0][i] for i in utrs5])
            end5 = max([exons[1][i] for i in utrs5])
        if utrs3:
            st3 = min([exons[0][i] for i in utrs3])
            end3 = max([exons[1][i] for i in utrs3])
        if utrs5 and utrs3:
            exons[0] = [st5] + [exons[0][i] for i in exonls] + [st3]
            exons[1] = [end5] + [exons[1][i] for i in exonls] + [end3]
            exons[2] = ['five_prime_UTR'] + [exons[2][i] for i in exonls] + ['three_prime_UTR']
        elif utrs5:
            exons[0] = [st5] + [exons[0][i] for i in exonls]
            exons[1] = [end5] + [exons[1][i] for i in exonls]
            exons[2] = ['five_prime_UTR'] + [exons[2][i] for i in exonls]
        elif utrs3:
            exons[0] = [exons[0][i] for i in exonls] + [st3]
            exons[1] = [exons[1][i] for i in exonls] + [end3]
            exons[2] = [exons[2][i] for i in exonls] + ['three_prime_UTR']
        exons[2] = [x for _, x in sorted(zip(exons[0], exons[2]))]
        exons[1] = [x for _, x in sorted(zip(exons[0], exons[1]))]
        exons[0] = sorted(exons[0])
    return (genes, exon_pos)

def annotBPs(sv, bp_1, genes, exon_pos):
    if not bp_1[0] in list(genes.keys()) and not 'chr' + bp_1[0] in list(genes.keys()):
        sv.genes.append(())
    else:
        if bp_1[0] in list(genes.keys()):
            genels = genes[bp_1[0]]
        else:
            genels = genes['chr' + bp_1[0]]
        ind1 = bisect.bisect_right(genels[1], bp_1[1])
        ind2 = bisect.bisect_left(genels[2], bp_1[1])
        if ind1 == ind2+1:
            if genels[0][ind2] in list(exon_pos.keys()):
                exons = exon_pos[genels[0][ind2]]
                ind1a = bisect.bisect_right(exons[0], bp_1[1])
                ind1b = bisect.bisect_left(exons[1], bp_1[1])
                if ind1a == 0 or ind1b == len(exons):
                    ty = 'promoter/utr'
                elif ind1a == ind1b:
                    ty = 'intron' + exons[2][ind1a - 1][-1]
                else:
                    ty = exons[2][ind1a - 1]
            sv.genes.append((genels[0][ind2], ty, exons[3][0]))
        else:
            sv.genes.append(())
            
def annot_SVS(genes, exon_pos, svs, by_gene):
    for sv in svs:
        annotBPs(sv, sv.bp_1, genes, exon_pos)
        annotBPs(sv, sv.bp_2, genes, exon_pos)
    for sv in svs:
        if not sv.genes[1] and not sv.genes[0]:
            sv.genes.append('')
            continue
        if not sv.genes[1] or not sv.genes[0]:
            sv.genes.append('Gene-NonCoding')
            continue
        if sv.genes[0] == sv.genes[1]:
            by_gene[sv.genes[0][0]].append(sv)
            if 'exon'  in sv.genes[0][1]:
                svlen = sv.bp_2[1] - sv.bp_1[1]
                if svlen // 3:
                    sv.genes.append('frameshift')
                else:
                    sv.genes.append('withinexon')
            elif 'intron' in sv.genes[0][1]:
                sv.genes.append('intronic')
            else:
                sv.genes.append('promoter/UTR')
        else:
            if sv.genes[0][0] == sv.genes[1][0]:
                by_gene[sv.genes[0][0]].append(sv)
                if not sv.genes[0][1] == sv.genes[1][1]:
                    sv.genes.append('between_exons')
            elif not sv.genes[0][0] == sv.genes[1][0]:
                by_gene[sv.genes[0][0]].append(sv)
                by_gene[sv.genes[1][0]].append(sv)
                if sv.direction_1 == sv.direction_2 and not sv.genes[0][2] == sv.genes[1][2]:
                    sv.genes.append('possible_fusion')
                elif not sv.direction_1 == sv.direction_2 and sv.genes[0][2] == sv.genes[1][2]:
                    sv.genes.append('possible_fusion')
                else:
                    sv.genes.append('between_genes')
                    
def annot_CNAs(genes, cnas, ploidy, by_gene):
    hps = [1,2]
    for ref_id, genels in genes.items():
        if not (ref_id,1) in cnas.keys():
            continue
        for hp in hps:
            cn_prof = cnas[ref_id, hp]
            cn_start = [cn.pos_1 for cn in cn_prof]
            cn_end = [cn.pos_2 for cn in cn_prof]
            cn = [cn.cn for cn in cn_prof]
            for i, gene in enumerate(genels[0]):
                ind1 = bisect.bisect_right(cn_start, genels[1][i])
                ind2 = bisect.bisect_left(cn_end, genels[2][i])
                by_gene[gene].append(ref_id + ':' + str(genels[1][i]) + '-' + str(genels[2][i]))
                if not ind1 - ind2 == 1:
                    cn_prof[ind1-1].genes.append((gene, 'disturbed'))
                    by_gene[gene].append(cn_prof[ind1-1])
                elif cn[ind1-1] > ploidy[hp-1]:
                    cn_prof[ind1-1].genes.append(gene)
                    cn_prof[ind1-1].dir1 = 'AMP'
                    by_gene[gene].append(cn_prof[ind1-1])
                elif cn[ind1-1] < ploidy[hp-1]:
                    cn_prof[ind1-1].genes.append(gene)
                    cn_prof[ind1-1].dir1 = 'DEL'
                    by_gene[gene].append(cn_prof[ind1-1])

def check_complexSV(cnas, svs):
    svls = defaultdict(list)
    for sv in svs:
         svls[sv.vcf_id] = sv

def run_command(cmd):
    p = subprocess.Popen(cmd, shell= True)
    p.wait()

def write_ins(svs, ref, t, specie):
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
    run_command(f"RepeatMasker -species {specie} temp_ins.fa")
    run_command(f"minimap2 -ax map-ont {ref} temp_ins.fa -k 17 -y -K 5G -t {t} --eqx | samtools sort -@ {t} -m 4G > temp_ins.bam")
    run_command(f"samtools index -@ {t} temp_ins.bam")
    

def get_repeat(svls):
    rep_file = open('temp_ins.fa.out', 'r')
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
    THR = 25
    bpls = []
    f = open("temp_bpseq.tsv", "w")
    for sv in svs:
        if sv.sv_type == 'INS':
            continue
        bpls.append(sv)
        run_command(f'samtools faidx {ref} {sv.bp_1[0]}:{max(0, sv.bp_1[1]-THR)}-{sv.bp_1[1]+THR} >> temp_bpseq.tsv')
        run_command(f'samtools faidx {ref} {sv.bp_2[0]}:{max(0, sv.bp_2[1]-THR)}-{sv.bp_2[1]+THR} >> temp_bpseq.tsv')
    
    f.close()
  
    f = open('temp_bpseq.tsv', 'r')
    t = 0
    seq = []
    tt = 0
    for line in f:
        if line.startswith('>'):
            continue
        else:
            seq.append(line.strip())
            t+=1
        if t == 2:
            check_homology(seq, bpls[tt])
            t = 0
            seq = []
            tt+=1

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
        
def annot_ins(svs, ref,t, rm_bed, specie):
    svls = defaultdict(list)
    for sv in svs:
         svls[sv.vcf_id] = sv
    write_ins(svs, ref, t, specie)
    get_repeat(svls)
    get_align(svls)
    get_tel(svs)
    annot_bp_repeat(svls, rm_bed)

def output_svs(svs, out_dir):
    fopen = open(out_dir + '/annotated_svs.tsv', 'w')
    header = '\t'.join(['chr_1', 'pos_1', 'chr_2', 'pos_2', 'n_supp', 'vaf', 'strand',\
                         'gene_name1', 'pos_1', 'strand_1', 'Repeat_Anno1',  'gene_name2', 'pos_2', \
                             'strand_2', 'Repeat_Anno2', 'Type', 'Microhomology', 'VNTR', \
                                 'Telomere_repeat', 'Aligned_pos(INS)', 'repeat_annot'])
    fopen.write(header)
    fopen.write('\n')
    for sv in svs:
        fopen.write(sv.to_str())
        fopen.write('\n')
    fopen.close()
    
def output_genes(by_gene, out_dir):
    f = open(out_dir + '/by_gene.tsv', 'w')
    header = '\t'.join(['Gene Name','Gene Pos', 'HP1', 'HP1 CN', 'HP2', 'HP2 CN', 'SV'])
    f.write(header)
    f.write('\n')
    for gene, val in by_gene.items():
        cnals = [v for v in val if type(v) is CNA]
        svls = [v for v in val if type(v) is SV]
        pos = [v for v in val if type(v is str)][0]
        if not svls and not cnals:
            continue
        cna_ls = [' \t ', ' \t ']
        sv_ls = []
        for cna in cnals:
            cna_ls[cna.haplotype-1] = cna.dir1 + '\t' + str(cna.cn)
        for sv in svls:
            sv_ls.append(sv.genes[-1])
        f.write(gene + '\t' + pos + '\t'+ '\t'.join(cna_ls) + '\t' + ';'.join(sv_ls))
        f.write('\n')
    f.close()


def annotate_things(args):
    cna_vcf, vcf_file, t, ref, out_dir = args.cna_vcf, args.vcf_file, args.threads, args.ref, args.out_dir
    by_gene = defaultdict(list)
    (genes, exon_pos) = get_genes(args.gff_file)
    svs = get_bps(vcf_file)
    cnas = []
    if cna_vcf:
        cnas, ploidy = get_CNA(cna_vcf, svs)
        annot_CNAs(genes, cnas, ploidy, by_gene)
    annot_SVS(genes, exon_pos, svs, by_gene)
    annot_ins(svs, ref,t, args.rm_bed, args.specie)
    get_microhomology(svs, ref)
    output_svs(svs, out_dir)
    output_genes(by_gene, out_dir)
    return (genes, cnas, exon_pos, svs, by_gene)
    

