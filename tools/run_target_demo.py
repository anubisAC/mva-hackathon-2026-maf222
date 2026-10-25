"""Reproduce the unfitted virtual target grid with no patient inputs or network."""
import argparse,csv,gzip,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from mva.target_first import sample_parameters, simulate
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,default=Path('outputs/public_target_demo'))
a=p.parse_args()
if a.output.exists(): raise SystemExit('Choose a fresh output directory.')
a.output.mkdir(parents=True)
policy=json.loads(Path('config/target_first_framework.json').read_text())
sim=policy['simulation']
samples=sample_parameters(sim)
contexts=[('healthy_tissue_rescue','same_residual_state')]+[('tumor_control',s) for s in sim['tumor_states']]
count=0
with gzip.open(a.output/'virtual_experiments.csv.gz','wt',encoding='utf-8',newline='') as stream:
 writer=None
 for parameters in samples:
  for arm,state in contexts:
   baseline=simulate(parameters,'vehicle',0,sim['time'],state)
   for intervention in sim['interventions']:
    for strength in sim['strengths']:
     result=simulate(parameters,intervention,strength,sim['time'],state)
     row={'scenario':parameters['scenario'],'sample':parameters['sample'],'arm':arm,'tumor_state':state if arm=='tumor_control' else '', 'intervention':intervention,'strength':strength,**result}
     for metric,value in result.items(): row['delta_'+metric]=value-baseline[metric] if value is not None and baseline[metric] is not None else None
     if writer is None: writer=csv.DictWriter(stream,fieldnames=list(row));writer.writeheader()
     writer.writerow(row);count+=1
(a.output/'summary.json').write_text(json.dumps({'evaluations':count,'measured_drug_profiles':0,'fitted':False,'interpretation':'Analyst-selected molecular sensitivity scenarios; not patient, tumor-survival or drug-response predictions.'},indent=2))
print(f'Completed {count:,} unfitted virtual evaluations; zero measured drug profiles.')
