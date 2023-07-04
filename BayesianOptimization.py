# this file iteratively performs Bayesian optimization of a function to find the global minima
# the following functions exist here:
# ->
#
# source code from:
# https://josh-melnick.medium.com/1d-and-2d-gpyopt-bayesian-optimization-of-real-world-experiments-explicit-example-code-fc4a377b48ff
# ^ truly a blessing good sir

# imports:
import numpy as np
import GPyOpt


# functions:
def next_node(datapath):
    data = read_in_bay_opt_param(datapath)
    inputs, outputs = format_data(data)

    xy_init = inputs
    z_init = outputs # * (-1)

    domain = [{'name': 'sr', 'type': 'continuous', 'domain': (0.8, 1.1)},
              {'name': 'd', 'type': 'continuous', 'domain': (10, 70)}]

    bo_step = GPyOpt.methods.BayesianOptimization(
        f=None,
        domain=domain,
        model_type='GP',                # Gaussian
        acquisition_type='MPI',       # Maximum Probability of Improvement
        # acquisition_type='EI',          # Expected Improvement
        X=xy_init,
        Y=z_init
        )

    x_next = bo_step.suggest_next_locations()

    print("Value of (x,y) that next should be evaluated, to minimize the objective func:"+str(x_next))
    bo_step.plot_acquisition()


def read_in_bay_opt_param(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
    ret = []
    for ln in lines:
        temp = [float(x) for x in ln.split()]
        ret.append([temp[0], temp[1], temp[2]])
    return ret


def format_data(data):
    inputs = np.array([[0, 0]])
    outputs = np.array([[0]])
    for d in data:
        inputs = np.append(inputs, np.array([[d[0], d[1]]]), axis=0)
        outputs = np.append(outputs, np.array([[d[2]]]), axis=0)
    inputs = np.delete(inputs, 0, 0)
    outputs = np.delete(outputs, 0, 0)
    return inputs, outputs
