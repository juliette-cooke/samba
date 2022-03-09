"""
Functions for preparing reactions in a metabolic model.
"""


def set_exchanges_rxn_bounds(model, exchange_min):
    """

    :param model:
    :param exchange_min:
    :return: model:
    """
    # Set the exchange reaction bounds
    for rxn in model.exchanges:
        rxn.lower_bound = -exchange_min
        rxn.upper_bound = 1000
    return model


def parse_ko_rxns(ko_filename):
    with open(ko_filename, 'r') as ko_file:
        ids_to_knockout = [i.strip().split(' ') for i in ko_file]
    return ids_to_knockout
