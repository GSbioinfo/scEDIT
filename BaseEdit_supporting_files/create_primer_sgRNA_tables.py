import pandas as pd 
from pyfaidx import Fasta

Tap_primers_mapping= pd.read_csv("blast_tapestri_primers.bsl",delimiter="\t",header=None)
