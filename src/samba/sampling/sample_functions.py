"""Functions related to sampling."""
import time
from cobra.sampling import sample
import logging
log = logging.getLogger(__name__)


def sample_time(model, nsamples, processors):
    start_time = time.time()
    s = sample(model, nsamples, processes=processors, )
    elapsed_time = time.time() - start_time
    log.info("Total elapsed time: {:.2f} sec".format(elapsed_time))
    return s
