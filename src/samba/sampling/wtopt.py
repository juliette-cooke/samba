"""Functions for optimising the WT."""


def wtopt_r2(rxn, model, model_clean, eps):
    with model_clean:  # Keeps the modifications to within the with
        # change the objective to the reaction:
        model_clean.objective = rxn
        # optimize for that objective
        FBA = model_clean.optimize(objective_sense="maximize")
        FBA_sol = FBA.objective_value
        print("FBA sol: " + str(FBA_sol))
        # set the new bounds as a percentage of the solution
        if FBA_sol > 0:
            rxn.lower_bound = FBA_sol * eps
        else:
            with model_clean:
                model_clean.objective = rxn
                FBA = model_clean.optimize(objective_sense="minimize")
                FBA_sol = FBA.objective_value
            # print("Min: " + str(FBA_sol))
            rxn.lower_bound = FBA_sol * eps
            rxn.upper_bound = 1000  # redundant? questionable? what about backwards reactions?
    return rxn.lower_bound, rxn.upper_bound
