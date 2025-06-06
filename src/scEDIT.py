import argparse
import BaseEdit
import InDel
from pathlib import Path

def main():
    """
    Usage: python scEDIT.py --mode BaseEdit/InDel --samples path/SampleList.xlsx --outDir path/output_directory

    Single cell Edit detection and Identification Tool (scEDIT) runs two modes of analysis: InDel and BaseEdit
        -InDel mode is used for single cell analysis of CRISPR-Cas induced indels 
        -BaseEdit mode is used for single cell analysis of CRISPR-Cas induced base edits
    
    Each mode required different set of excel files containing amplicon primer information and gRNA details.
    
    The excel files need to have mode specific format.
    
    For more detail refer to the example files provide in github repo: 
  
    """
    parser = argparse.ArgumentParser(description="Run selected module")
    parser.add_argument("--mode", choices=["BaseEdit", "InDel"], required=True, help="Which mode to run")
    parser.add_argument("--samples", required=True, type=str, help="Sample sheet name with full path")
    parser.add_argument("--outDir", required=True, type=str, help="output directory full path")
    args = parser.parse_args()
    scriptdir = str(Path(__file__).resolve().parent)
    print(f"Script is located in: {scriptdir}")
    bindir = scriptdir+"/bin"
    inout_dict = {
    "samples": args.samples,
    "outDir": args.outDir,
    "scriptdir": scriptdir,
    "binDir": bindir
    }
    if args.mode == "BaseEdit":
        BaseEdit.run(inout_dict)
        print("selected BaseEdit mode")
    elif args.mode == "InDel":
        print("selected InDel mode")
        InDel.run(inout_dict)
        
if __name__ == "__main__":
    main()
