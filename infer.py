import numpyro
import numpyro.distributions as dist
import jax.numpy as jnp
from numpyro.infer import MCMC, NUTS
from numpyro.infer import Predictive
from enum import Enum

class kSeroStatus(Enum):
    NEGATIVE = 0
    POSITIVE = 1
    UNKNOWN = 2

class kSeenStatus(Enum):
    UNSEEN = 0
    SEEN = 1
    CAPTURED = 2

#Holds time series data collected for a capture mark recapture bird
class TaggedBird():
    time_tagged = None
    sero_status = jnp.array([])
    seen_status = jnp.array([])
    island_number = jnp.array([])
    def __init__(sero_series, cmr_series, island_series, start_time):
        sero_status = sero_series
        seen_status = cmr_series
        island_number = island_series
        time_tagged = start_time

#generative model that will also be used for inference
def generate_cmr_data_set(num_sites):
    numpyro.sample