"""Run the complete project without opening plotting windows."""
from pathlib import Path
import argparse
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from experiments import run_all

if __name__=='__main__':
    parser=argparse.ArgumentParser(description='Coherent RTD transport and numerical validation')
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'figures')
    args=parser.parse_args()
    result=run_all(args.output)
    print('Saved figures, data and summary to',args.output.resolve())
    print('First peak/valley:',result['baseline_5mV'])
    print('Checks:',result['checks'])
