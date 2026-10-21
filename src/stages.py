"""Small experiments for learning, without running the entire study."""
from pathlib import Path
import argparse
from experiments import baseline, sensitivity, verification, sampling_failure
from negf import Device

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('stage',choices=['baseline','verify','sensitivity','sampling'])
    p.add_argument('--output',type=Path,default=Path('reproduction'))
    a=p.parse_args(); a.output.mkdir(parents=True,exist_ok=True)
    functions={'baseline':baseline,'verify':verification,'sensitivity':sensitivity,'sampling':sampling_failure}
    print(functions[a.stage](a.output,Device()))
