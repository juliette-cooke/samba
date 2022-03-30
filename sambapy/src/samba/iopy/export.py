"""Functions for exporting and writing files."""
import os
from setup.prepare_reactions import parse_rxns


def write_sampling(s_results, out_path, model_name, n_samples, type):
    final_path = str(out_path) + str(os.path.basename(os.path.splitext(model_name)[0])) + "_sampling_" + str(
                n_samples) + "_" + type + ".csv.gz"
    s_results.to_csv(final_path, index=False, compression='gzip')
    print("Wrote to "+final_path)


def extract_results(s, results,):
    # Extract from results
    if results is None:
        results_columns = [col for col in s if col.startswith('EX_')]
    else:
        results_columns = parse_rxns(results)
    s_results = s[results_columns].round(3)
    return s_results


def export_metab_dict(model):
    # Create a metabolite ID to name dict for plotting in R
    metab_dict = {}
    for ex_rxn in model.exchanges:
        if len(next(iter(ex_rxn.metabolites)).name) < 30:
            metab_dict[ex_rxn.id] = next(iter(ex_rxn.metabolites)).name
        else:
            metab_dict[ex_rxn.id] = next(iter(ex_rxn.metabolites)).id[:-3]
    return metab_dict
