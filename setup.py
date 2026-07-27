from setuptools import setup


SCRIPT_FILES = [
    "scripts/FPKM-UQ.py",
    "scripts/FPKM_count.py",
    "scripts/RNA_fragment_size.py",
    "scripts/RPKM_saturation.py",
    "scripts/bam2fq.py",
    "scripts/bam2wig.py",
    "scripts/bam_stat.py",
    "scripts/clipping_profile.py",
    "scripts/deletion_profile.py",
    "scripts/divide_bam.py",
    "scripts/geneBody_coverage.py",
    "scripts/geneBody_coverage2.py",
    "scripts/infer_experiment.py",
    "scripts/inner_distance.py",
    "scripts/insertion_profile.py",
    "scripts/junction_annotation.py",
    "scripts/junction_saturation.py",
    "scripts/mismatch_profile.py",
    "scripts/normalize_bigwig.py",
    "scripts/overlay_bigwig.py",
    "scripts/read_GC.py",
    "scripts/read_NVC.py",
    "scripts/read_distribution.py",
    "scripts/read_duplication.py",
    "scripts/read_hexamer.py",
    "scripts/read_quality.py",
    "scripts/sc_bamStat.py",
    "scripts/sc_editMatrix.py",
    "scripts/sc_seqLogo.py",
    "scripts/sc_seqQual.py",
    "scripts/split_bam.py",
    "scripts/split_paired_bam.py",
    "scripts/tin.py",
]


if __name__ == "__main__":
    setup(scripts=SCRIPT_FILES)
