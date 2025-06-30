#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 14:05:48 2023

@author: keskusa2
"""

#!/usr/bin/env python3

import sys
import shutil
import pysam
import argparse
import os
from multiprocessing import Pool
from collections import defaultdict,Counter
import logging

from padfoot.annot import annotate_things
from padfoot.preprocess import generate_gff, generate_rm
from padfoot_cluster_cns import cluster_cn
from padfoot.__version__ import __version__


logger = logging.getLogger()


def _enable_logging(log_file, debug, overwrite):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    if not debug:
        console_log.setLevel(logging.INFO)

    if overwrite:
        open(log_file, "w").close()
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)

def _version():
    return __version__


def main():
    
    SAMTOOLS_BIN = "samtools"

    parser = argparse.ArgumentParser \
        (description="Find breakpoints and build breakpoint graph from a bam file")

    parser.add_argument("-v", "--version", action="version", version=_version())
    parser.add_argument("--severus_svs", dest="vcf_file",
                        metavar="path", required=True, default=None,
                        help="path to severus vcf file")
    parser.add_argument("--HP1", dest="hp1_file",
                        metavar="path", required=True, default=None,
                        help="path to Wakhan HP1 file")
    parser.add_argument("--HP2", dest="hp2_file",
                        metavar="path", required=True, default=None,
                        help="path to Wakhan HP2 file")
    parser.add_argument("--LOH", dest="loh_file",
                        metavar="path", required=True, default=None,
                        help="path to Wakhan LOH file")
    parser.add_argument("--ref", "-r", dest="ref",
                        metavar="path", required=True, default=None,
                        help="path to reference file")
    parser.add_argument("-t", "--threads", dest="threads",
                        default=8, metavar="int", type=int, help="number of parallel threads [8]")
    parser.add_argument("--out-dir", dest="out_dir",
                        default=None, required=True,
                        metavar="path", help="Output directory")
    parser.add_argument("--gff", dest="new_gff",
                        default=None,
                        metavar="path", help="GFF annotation file")
    parser.add_argument("--genome", dest="genome",
                        default='hg38', metavar="int", help="Either hg38, chm13, mm10")
    parser.add_argument("--rm", dest="rm_file",
                        default = None,
                        metavar="path", help="Repeat masker annotation file")
    parser.add_argument("-t", "--threads", dest="threads",
                        default=8, metavar="int", type=int, help="number of parallel threads [8]")
    parser.add_argument("--specie", dest="specie",
                        default=None, help="Specie")
    args = parser.parse_args()
    

    temp_dir = args.out_dir+ '/temp'
    if not os.path.isdir(args.out_dir):
        os.makedirs(temp_dir)
        
    if not  os.path.isdir(temp_dir):
        os.makedirs(temp_dir)
    
    os.chdir(args.out_dir+ '/temp')
    log_file = os.path.join(args.out_dir, "padfoot.log")
    _enable_logging(log_file, debug=False, overwrite=True)

    logger.info("Starting Padfoot " + _version())
    logger.debug("Cmd: %s", " ".join(sys.argv))
    logger.debug("Python version: " + sys.version)

    if args.new_gff:
        gff_path = args.out_dir+ 'gff_file.gff3'
        logger.debug("Generating new gff file: " + gff_path)
        generate_gff(args.new_gff, gff_path)
        args.gff_file = gff_path
    
    if args.rm_file:
        rm_path = args.out_dir+ 'rm.bed'
        logger.debug("Generating new rm file: " + rm_path)
        generate_rm(args.rm_file, rm_path)
        args.rm_file = rm_path
        
    if not args.rm_file and not args.genome:
        logger.error('Error: Please provide a repeat masker bed file or select a genome [hg38, mm10, chm13]')
    if not args.rm_file and args.genome:
        args.gff_file = 'beds/' + args.genome +'_rm.bed'
    
    if not args.new_gff and not args.genome:
        logger.error('Error: Please provide a gff file or select a genome [hg38, mm10, chm13]')
    if not args.new_gff and args.genome:
        args.gff_file = 'beds/' + args.genome +'.gff3.gz'
    
    if not args.genome and not args.specie:
        logger.error('Error: Please provide specie for repeat annotation')
        
    args.specie = 'mouse' if args.genome == 'mm10' else 'human'
    
    (genes, cnas, exon_pos, svs, by_gene) = annotate_things(args)
    cluster_cn(cnas, args.out_dir)
    

