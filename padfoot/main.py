#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 14:05:48 2023

@author: keskusa2
"""

#!/usr/bin/env python3

import sys
import shutil
import argparse
import os
import importlib
import logging

from padfoot.annot import annotate_things
from padfoot.preprocess import generate_rm
from padfoot.annotation import build_annotation, is_compact_annotation
from padfoot.copynumber import SEXES, GAIN_HIGH_MULTIPLIER
#from padfoot_cluster_cns import cluster_cn
from padfoot.__version__ import __version__


logger = logging.getLogger()

def _version():
    return __version__

def _check_python_dependencies():
    required_modules = {
        "pysam": "pysam",
        "pandas": "pandas",
        "numpy": "numpy",
        "Bio": "biopython",
    }

    missing = []
    for module_name, package_name in required_modules.items():
        try:
            importlib.import_module(module_name)
        except ImportError:
            missing.append(package_name)

    if missing:
        logger.error(
            "Missing required Python packages: %s",
            ", ".join(sorted(set(missing)))
        )
        logger.error(
            "Install them with conda, for example:\n"
            "  conda env create -f environment.yml\n"
            "or\n"
            "  conda install -c bioconda -c conda-forge %s",
            " ".join(sorted(set(missing)))
        )
        sys.exit(1)



def _check_external_dependencies(args, ref_path=None):
    required_bins = ["samtools", "minimap2"]
    missing_bins = [tool for tool in required_bins if shutil.which(tool) is None]

    if missing_bins:
        logger.error(
            "Missing required external tools in PATH: %s",
            ", ".join(missing_bins)
        )
        logger.error(
            "Please install them in your conda environment, e.g.:\n"
            "  conda install -c bioconda -c conda-forge %s",
            " ".join(missing_bins)
        )
        sys.exit(1)
    
    required_bins = ["RepeatMasker"]
    missing_bins = [tool for tool in required_bins if shutil.which(tool) is None]
    if missing_bins:
        args.run_repeatmasker = False
        logger.info('RepeatMasker is not found... \n Skipping RepeatMasker')

    if ref_path is not None:
        if not os.path.exists(ref_path):
            logger.error("Reference FASTA not found: %s", ref_path)
            sys.exit(1)

        fai_path = ref_path + ".fai"
        if not os.path.exists(fai_path):
            logger.error(
                "Reference FASTA index not found: %s\n"
                "Create it with:\n"
                "  samtools faidx %s",
                fai_path,
                ref_path,
            )
            sys.exit(1)
            
def _enable_logging(log_file, debug, overwrite):
    log_formatter = logging.Formatter(
        "[%(asctime)s] %(name)s: %(levelname)s: %(message)s",
        "%Y-%m-%d %H:%M:%S"
    )
    console_formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)s: %(message)s",
        "%Y-%m-%d %H:%M:%S"
    )

    logger.handlers.clear()

    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    console_log.setLevel(logging.DEBUG if debug else logging.INFO)

    if overwrite:
        open(log_file, "w").close()

    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler) 
    
    
def main():
    parser = argparse.ArgumentParser(
        description="Annotate structural variants and copy-number alterations"
    )
    parser.add_argument("-v", "--version", action="version", version=_version())
    parser.add_argument("--severus-vcf", dest="vcf_file", metavar="path", required=True,
                        help="path to Severus vcf file")
    parser.add_argument("--wakhan-profile", dest="wakhan_profile", metavar="path", required=True,
                        help="Wakhan (>= 0.5.0) integer_profile.bed of the chosen solution")
    parser.add_argument("--sex", dest="sex", choices=SEXES, default="unknown",
                        help="sex of the sample; no LOH is called on chrX of male samples [unknown]")
    parser.add_argument("--gain-high-multiplier", dest="gain_high_multiplier", default=GAIN_HIGH_MULTIPLIER,
                        type=float, metavar="float",
                        help="GAIN_HIGH when the gene's total CN exceeds this multiple of the modal CN [3]")
    parser.add_argument("--ref", "-r", dest="ref", metavar="path", required=True,
                        help="path to reference file")
    parser.add_argument("--out-dir", dest="out_dir", required=True, metavar="path",
                        help="Output directory")
    parser.add_argument("--gff", dest="new_gff", default=None, metavar="path",
                        help="gene annotation: a GFF3 (GENCODE/Ensembl style) or a table made by "
                             "scripts/build_annotation.py [bundled GENCODE annotation for --genome]")
    parser.add_argument("--genome", dest="genome", default="hg38", metavar="name",
                        help="Either hg38, chm13, mm10")
    parser.add_argument("--rm", dest="rm_file", default=None, metavar="path",
                        help="Repeat masker annotation file")
    parser.add_argument("-t", "--threads", dest="threads", default=8, type=int, metavar="int",
                        help="number of parallel threads [8]")
    parser.add_argument("--specie", dest="specie", default="human", help="Specie")
    parser.add_argument("--cancer-genes", dest="cancer_genes", default=None, metavar="path",
                        help="cancer gene list (TSV with Gene_symbol, Role_in_cancer, fusion columns) "
                             "[bundled list for --specie human, none otherwise]")
    parser.add_argument("--whole-gene-max-sv", dest="whole_gene_max_sv", default=1_000_000, type=int, metavar="int",
                        help="label genes lying entirely inside a DEL/DUP/INV up to this length [1000000]")
    parser.add_argument("--skip_RepeatMasker", dest="run_repeatmasker", action = "store_false", help="Skip RepeatMasker [True]")
    args = parser.parse_args()

    beds = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "beds"))

    args.out_dir = os.path.abspath(args.out_dir)
    args.vcf_file = os.path.abspath(args.vcf_file)
    args.wakhan_profile = os.path.abspath(args.wakhan_profile)
    args.ref = os.path.abspath(args.ref)

    os.makedirs(args.out_dir, exist_ok=True)

    temp_dir = os.path.join(args.out_dir, "temp")
    os.makedirs(temp_dir, exist_ok=True)
    os.chdir(temp_dir)

    log_file = os.path.join(args.out_dir, "padfoot.log")
    _enable_logging(log_file, debug=True, overwrite=True)
    _check_python_dependencies()
    _check_external_dependencies(args, args.ref)

    logger.info("Starting Padfoot %s", _version())
    logger.debug("Command: %s", " ".join(sys.argv))
    logger.debug("Python version: %s", sys.version)
    logger.info("Working directory: %s", os.getcwd())
    logger.info("Temporary directory: %s", temp_dir)
    logger.info("Output directory: %s", args.out_dir)

    logger.info("Input Severus VCF: %s", args.vcf_file)
    logger.info("Input Wakhan profile: %s", args.wakhan_profile)
    logger.info("Sex: %s", args.sex)
    logger.info("Reference FASTA: %s", args.ref)
    logger.info("Threads: %s", args.threads)
    logger.info("Requested species: %s", args.specie)
    logger.info("Requested genome preset: %s", args.genome)

    if not os.path.exists(args.vcf_file):
        logger.error("Severus VCF does not exist: %s", args.vcf_file)
        sys.exit(1)

    if not os.path.isfile(args.wakhan_profile):
        logger.error("Wakhan profile does not exist or is not a file: %s", args.wakhan_profile)
        sys.exit(1)
    if args.sex == "unknown":
        logger.warning("--sex not given: LOH is called on chrX as for a female sample")

    if not os.path.exists(args.ref):
        logger.error("Reference FASTA does not exist: %s", args.ref)
        sys.exit(1)

    if not os.path.exists(args.ref + ".fai"):
        logger.error("Reference FASTA index is missing: %s.fai", args.ref)
        logger.error("Create it with: samtools faidx %s", args.ref)
        sys.exit(1)

    logger.info("Reference FASTA index found: %s.fai", args.ref)

    if args.new_gff:
        if not os.path.exists(args.new_gff):
            logger.error("Annotation file does not exist: %s", args.new_gff)
            sys.exit(1)
        if is_compact_annotation(args.new_gff):
            args.gff_file = os.path.abspath(args.new_gff)
            logger.info("Using prebuilt annotation table: %s", args.gff_file)
        else:
            args.gff_file = os.path.join(args.out_dir, "gene_annotation.tsv.gz")
            logger.info("Converting custom GFF3 %s to %s", args.new_gff, args.gff_file)
            counts = build_annotation(args.new_gff, args.gff_file)
            logger.info("Annotation table: %d genes, %d transcripts", counts["genes"], counts["transcripts"])
    else:
        if not args.genome:
            logger.error("Please provide --gff or select --genome [hg38, mm10]")
            sys.exit(1)
        args.gff_file = os.path.join(beds, f"{args.genome}_genes.tsv.gz")
        if not os.path.exists(args.gff_file):
            logger.error("No bundled gene annotation for genome %s; provide one with --gff", args.genome)
            sys.exit(1)
        logger.info("Using bundled gene annotation: %s", args.gff_file)

    if args.rm_file:
        rm_path = os.path.join(args.out_dir, "rm.bed")
        logger.info("Using custom RepeatMasker input: %s", args.rm_file)
        logger.debug("Generating normalized RepeatMasker BED: %s", rm_path)
        generate_rm(args.rm_file, rm_path)
        args.rm_file = rm_path
        logger.info("Prepared RepeatMasker annotations: %s", args.rm_file)
    else:
        if not args.genome:
            logger.error("Please provide --rm or select --genome [hg38, mm10, chm13]")
            sys.exit(1)
        args.rm_file = os.path.join(beds, f"{args.genome}_rm.bed.gz")
        logger.info("Using bundled RepeatMasker annotations: %s", args.rm_file)

    if not args.specie:
        logger.error("Please provide --specie for repeat annotation")
        sys.exit(1)

    if args.genome == "mm10" and args.specie != "mouse":
        logger.warning("Genome preset mm10 usually implies species 'mouse'; current value: %s", args.specie)

    if args.cancer_genes:
        if not os.path.exists(args.cancer_genes):
            logger.error("Cancer gene list does not exist: %s", args.cancer_genes)
            sys.exit(1)
        args.cancer_genes = os.path.abspath(args.cancer_genes)
    elif args.specie == "human":
        args.cancer_genes = os.path.join(beds, "cancer_genes.tsv")
    if args.cancer_genes:
        logger.info("Cancer gene list: %s", args.cancer_genes)
    else:
        logger.info("No cancer gene list for species %s; cancer roles and known fusions not annotated", args.specie)

    logger.info("Annotation inputs resolved successfully")
    logger.info("Gene annotation: %s", args.gff_file)
    logger.info("RepeatMasker file: %s", args.rm_file)

    logger.info("Starting annotation pipeline")
    annotation, profile, svs, by_gene = annotate_things(args)

    logger.info("Annotation pipeline finished")
    logger.info("Summary: genes=%d, cn_fragments=%d, svs=%d, by_gene=%d",
                len(annotation.genes), len(profile) if profile else 0, len(svs), len(by_gene))
    logger.info("Log file written to: %s", log_file)
    logger.info("Padfoot completed successfully")