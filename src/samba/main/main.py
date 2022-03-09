from json import dump
from sambapy.src.samba.io.read_model import import_model
from sambapy.src.samba.setup.prepare_reactions import set_exchanges_rxn_bounds, parse_ko_rxns
from sambapy.src.samba.sampling.wtopt import wtopt_r2
import logging
import argparse
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
    parser.add_argument("-q", "--quiet", action="store_true", help="Use this flag to silence INFO logging.")
    parser.add_argument("--log", help="Log file path (without filename). Set to None to use the outfile path.")
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
    # Assign the INFO level to print
    print = logging.info
    if args.log is not None:
        logfile_path = args.log
    else:
        logfile_path = args.outpath
    # Set level to only warnings or higher if QUIET is true
    logging.basicConfig(level=logging.WARNING if QUIET else logging.INFO,
                        format='%(levelname)s:%(message)s',
                        filename=logfile_path + "out.log", encoding="utf-8")
    # For important warnings that will be printed regardless of QUIET
    iprint = logging.warning

    model_file = args.model

    # TEST
    model_file = "/home/juliette/these/data/models/test_samba/Human-GEM.xml"

    model = import_model(model_file)
    model.solver = args.solver

    eps = 0.05
    # TODO: make as parameter

    if args.exchangemin is not None:
        model = set_exchanges_rxn_bounds(model, args.exchangemin)

    # TODO: Check that all KO reactions exist in the model before running everything

    ids_to_knockout = parse_ko_rxns(args.ko)

    model_clean = model.copy()
    for rxn_group in ids_to_knockout:
        for r in rxn_group:
            rxn = model.reactions.get_by_id(r)
            rxn.lower_bound, rxn.upper_bound = wtopt_r2(rxn, model, model_clean)
        print(str(model.objective))
        fbag = model.optimize()  # The optimized objective value
        # The minimum objective value we want to reach:
        print("Obj value:" + str(fbag.objective_value))
        obj = args.biomass * fbag.objective_value  # Store this so we can use the same value for the KOs
        print(obj)
        model.reactions.biomass_reaction.lower_bound = obj

    # Tryin to set the existing objective coeff to 0: not sure if it works yet
    # obj_rxn = model.reactions.get_by_id(str(model.objective.expression.args[0])[4:])
    # obj = {obj_rxn: 0}
    # model.objective = obj