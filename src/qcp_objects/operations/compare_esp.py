from . import *

def load_maps(trial_map, reference_map):
    # Load the maps
    maps={
        'trial': trial_map,
        'reference':reference_map
    }
    for k,v in maps.items():
        try:
            maps[k]=esp_map_cl(v)
        except Exception as ex:
            raise Exception(f"Error for loading \'{k}\' map from file {v}: {ex}")
    return maps['trial'], maps['reference']

@validate_call
def compare_esp(trial_map:file, reference_map:file):
    
    trial_map, reference_map=load_maps(trial_map=trial_map, reference_map=reference_map)

    dif_map=trial_map.sub(reference_map)

    # could be generalized as method
    final_file=dif_map.to_file(onam='test.map.json',format='json')
    results=dict(
        files={
            'map':final_file
        },
        results={
            'statistical_indicators':dif_map.gen_statistical_indicators()
        }
    )
    # Drop the results
    results_file='results.yaml'
    with open(results_file, 'w') as wr:
        yaml.safe_dump(results,wr)
    print(f"Dropped results to {results_file}")
    return results_file


