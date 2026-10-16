"""Acquisition functions for the BoTorch sampler in Optuna (Appendix C).

Only EI_candidates_func (Expected Improvement) is used in the paper.
"""
from botorch.acquisition.analytic import ExpectedImprovement, ProbabilityOfImprovement, UpperConfidenceBound
from botorch.acquisition.monte_carlo import qExpectedImprovement
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from botorch.models.transforms.outcome import Standardize
from botorch.optim import optimize_acqf
from botorch.sampling import SobolQMCNormalSampler
from botorch.utils.transforms import normalize
from botorch.utils.transforms import unnormalize
from gpytorch.mlls import ExactMarginalLogLikelihood
import numpy as np
import optuna
import torch


def PI_candidates_func(
    train_x, train_obj, train_con, bounds, pending_x
):
        
    train_x = normalize(train_x, bounds=bounds)
    
    model = SingleTaskGP(train_x, train_obj, outcome_transform=Standardize(m=train_obj.size(-1)))
    mll = ExactMarginalLogLikelihood(model.likelihood, model)
    fit_gpytorch_mll(mll)

    acqf = ProbabilityOfImprovement(
        model=model,
        best_f=train_obj.max(),
    )

    standard_bounds = torch.zeros_like(bounds)
    standard_bounds[1] = 1

    candidates, _ = optimize_acqf(
        acq_function=acqf,
        bounds=standard_bounds,
        q=1,
        num_restarts=10,
        raw_samples=512,
        options={"batch_limit": 5, "maxiter": 200},
        sequential=True,
    )

    candidates = unnormalize(candidates.detach(), bounds=bounds)

    return candidates






def EI_candidates_func(
    train_x, train_obj, train_con, bounds, pending_x
):
        
    train_x = normalize(train_x, bounds=bounds)
    
    model = SingleTaskGP(train_x, train_obj, outcome_transform=Standardize(m=train_obj.size(-1)))
    mll = ExactMarginalLogLikelihood(model.likelihood, model)
    fit_gpytorch_mll(mll)

    acqf = ExpectedImprovement(
        model=model,
        best_f=train_obj.max(),
    )

    standard_bounds = torch.zeros_like(bounds)
    standard_bounds[1] = 1

    candidates, _ = optimize_acqf(
        acq_function=acqf,
        bounds=standard_bounds,
        q=1,
        num_restarts=10,
        raw_samples=512,
        options={"batch_limit": 5, "maxiter": 200},
        sequential=True,
    )

    candidates = unnormalize(candidates.detach(), bounds=bounds)

    return candidates




def LCB_candidates_func(
    train_x, train_obj, train_con, bounds, pending_x
):
        
    train_x = normalize(train_x, bounds=bounds)
    
    model = SingleTaskGP(train_x, train_obj, outcome_transform=Standardize(m=train_obj.size(-1)))
    mll = ExactMarginalLogLikelihood(model.likelihood, model)
    fit_gpytorch_mll(mll)

    acqf = UpperConfidenceBound( 
        model=model,
        beta=0.1,# width of the confidence interval
    )

    standard_bounds = torch.zeros_like(bounds)
    standard_bounds[1] = 1

    candidates, _ = optimize_acqf(
        acq_function=acqf,
        bounds=standard_bounds,
        q=1,
        num_restarts=10,
        raw_samples=512,
        options={"batch_limit": 5, "maxiter": 200},
        sequential=True,
    )

    candidates = unnormalize(candidates.detach(), bounds=bounds)

    return candidates




def LCB_logn_candidates_func(
    train_x, train_obj, train_con, bounds, pending_x
):
        
    train_x = normalize(train_x, bounds=bounds)
    
    model = SingleTaskGP(train_x, train_obj, outcome_transform=Standardize(m=train_obj.size(-1)))
    mll = ExactMarginalLogLikelihood(model.likelihood, model)
    fit_gpytorch_mll(mll)
    
    n_trials = train_x.shape[0] + 1 - 10 # n_trials includes the 10 initial trials
    beta = 0.1 * np.log(n_trials)
    
    acqf = UpperConfidenceBound(
        model=model,
        beta=beta,
    )

    standard_bounds = torch.zeros_like(bounds)
    standard_bounds[1] = 1

    candidates, _ = optimize_acqf(
        acq_function=acqf,
        bounds=standard_bounds,
        q=1,
        num_restarts=10,
        raw_samples=512,
        options={"batch_limit": 5, "maxiter": 200},
        sequential=True,
    )

    candidates = unnormalize(candidates.detach(), bounds=bounds)

    return candidates

