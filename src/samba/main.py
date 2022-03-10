import time
from json import dump
from iopy.read_model import import_model
from setup.prepare_reactions import set_exchanges_rxn_bounds, parse_rxns
from sampling.wtopt import wtopt_r2, optimise_wt
from sampling.sample_functions import sample_time
from iopy.export import write_sampling, extract_results
import logging
import argparse
from cobra.sampling import sample

MIN_VAL = 0
MAX_VAL = 1


def range_limited_float_type(arg):
    """ Type function for argparse - a float within some predefined bounds """
    try:
        f = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be a floating point number")
    if f < MIN_VAL or f > MAX_VAL:
        raise argparse.ArgumentTypeError("Argument must be < " + str(MAX_VAL) + "and > " + str(MIN_VAL))
    return f


if __name__ == '__main__':
    description = "Package for running flux sampling on metabolic models."

    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-m", "--model", help="Metabolic model. Supported file formats: SBML, json, mat.")
    parser.add_argument("-n", "--nsamples", type=int, help="The number of samples")
    parser.add_argument("-p", "--processors", type=int, help="Number of processors")
    parser.add_argument("-t", "--thinning", type=int, required=False, help="The thinning factor of the generated "
                                                                           "sampling chain. A thinning of 10 means "
                                                                           "samples are returned every 10 steps.")
    parser.add_argument("-o", "--outpath", help="Outfile path (without filename)")
    parser.add_argument("-k", "--ko", help="KO file containing reactions to KO, "
                                           "specify nothing if you want to sample WT", default=None)
    parser.add_argument("-r", "--results", help="File containing reactions to output, "
                                                "specify nothing if you want only exchange reactions", default=None)
    parser.add_argument("-q", "--quiet", action="store_true", help="Use this flag to silence INFO logging.")
    parser.add_argument("-d", "--debug", action="store_true", help="Use this flag to turn on DEBUG logging.")
    parser.add_argument("--dryrun", action="store_true", help="Use this flag to run the code without running sampling.")
    parser.add_argument("--log", help="Log file path + filename. Set to None to output to console.",
                        default=None)
    parser.add_argument("-b", "--biomass", type=range_limited_float_type,
                        help="Number between 0 and 1, fraction of biomass to optimize", default=0)
    parser.add_argument("--solver", help="Solver to use", choices=["cplex", "gurobi", "glpk"], default="cplex")
    parser.add_argument("--exchangemin", default=10,
                        help="The value used to set the minimum flux through exchange reactions (will be negative). "
                             "Set to None if you want the default exchange reaction bounds.")

    args = parser.parse_args()

    # Write parameter config file
    config_dict = args.__dict__
    with open(args.outpath + 'config.json', 'w') as config_file:
        dump(config_dict, config_file, indent=4)

    # Logging setup
    QUIET = args.quiet
    if args.log is not None:
        logfile_path = args.log
    else:
        logfile_path = None
    # Set level to only warnings or higher if QUIET is true
    if QUIET:
        level = logging.WARNING
    else:
        if args.debug:
            level = logging.DEBUG
        else:
            level = logging.INFO
    # Assign the INFO level to print
    print = logging.info
    logging.basicConfig(level=level,
                        format='%(levelname)s:%(message)s',
                        filename=logfile_path)
    # For important warnings that will be printed regardless of QUIET
    iprint = logging.warning

    dprint = logging.debug

    model_file = args.model

    # TEST
    # model_file = "/home/juliette/these/data/models/test_samba/RECON1.xml"

    model = import_model(model_file)
    model.solver = args.solver

    eps = 0.05
    # TODO: make as parameter

    if args.exchangemin is not None:
        model = set_exchanges_rxn_bounds(model, args.exchangemin)

    # TODO: Check that all KO reactions exist in the model before running everything

    ids_to_knockout = parse_rxns(args.ko)

    model_clean = model.copy()

    # If no KO ids were provided, just do a basic WT sampling
    if ids_to_knockout is None:
        if args.biomass != 0:
            model.reactions.biomass_reaction.lower_bound = optimise_wt(model, args.biomass)
        if not args.dryrun:
            s = sample_time(model, args.nsamples, args.processors)
            s_results = extract_results
            write_sampling(s_results, args.outpath, args.model, args.nsamples)
    else:
        # If KO ids were provided, loop over them and sample the WT then the KO
        for rxn_group in ids_to_knockout:
            # Set up the rxn bounds
            for r in rxn_group:
                rxn = model.reactions.get_by_id(r)
                rxn.lower_bound, rxn.upper_bound = wtopt_r2(rxn, model_clean, eps)
            model.reactions.biomass_reaction.lower_bound = optimise_wt(model, args.biomass)

            # SAMPLING WT
            # Now we need to run the WT sampling on the model with these bounds and objective value
            if not args.dryrun:
                s = sample_time(model, args.nsamples, args.processors)
                s_results = extract_results
                write_sampling(s_results, args.outpath, args.model, args.nsamples, "WT")

            with model:
                for r in rxn_group:
                    rxn = model.reactions.get_by_id(r)
                    # Set all the reactions to KO to 0
                    rxn.lower_bound = rxn.upper_bound = 0
                # Re-optimise for the biomass function now that the KO reaction bounds are 0
                model.reactions.biomass_reaction.lower_bound = optimise_wt(model, args.biomass)

                # SAMPLING KO
                # Now we need to run the KO sampling on the model with these bounds and objective value
                if not args.dryrun:
                    s = sample_time(model, args.nsamples, args.processors)
                    s_results = extract_results
                    write_sampling(s_results, args.outpath, args.model, args.nsamples, "KO")




    # Trying to set the existing objective coeff to 0: not sure if it works yet
    # obj_rxn = model.reactions.get_by_id(str(model.objective.expression.args[0])[4:])
    # obj = {obj_rxn: 0}
    # model.objective = obj
