from json import dump
import sambapy
from sambapy.src.samba.io.read_model import import_model
import logging
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter

if __name__ == '__main__':
    description = "Package for running flux sampling on metabolic models."

    parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)
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
                        filename=logfile_path+"out.log", encoding="utf-8")
    # For important warnings that will be printed regardless of QUIET
    iprint = logging.warning

    # TEST
    model_file = "/home/juliette/these/data/models/raw/RECON1.xml"
    model = import_model(model_file)
