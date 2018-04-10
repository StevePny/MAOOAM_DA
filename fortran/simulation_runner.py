"""
By Eviatar Bach, see https://github.com/eviatarbach/simulation_runner
"""
import subprocess
import itertools
import time
import math
import operator
from functools import reduce
from abc import ABC, abstractmethod

from mako.template import Template


class _Namer(ABC):
    def __init__(self, **kwargs):
        pass

    @abstractmethod
    def next(self, keys, param_set):
        pass


class Hasher(_Namer):
    def next(self, keys, param_set):
        return hash(''.join(map(str, param_set))) % 65535


class String(_Namer):
    def next(self, keys, param_set):
        return '_'.join('{key}={param}'.format(key=key, param=param) for key, param in zip(keys, param_set))


class Sequential(_Namer):
    def __init__(self, length):
        self.count = -1
        self.zfill = math.floor(math.log10(length - 1) + 1)

    def next(self, keys, param_set):
        self.count += 1
        return str(self.count).zfill(self.zfill)


def run_simulation(command, config_path, template_path=None,
                   template_text=None, single_parameters={},
                   sweep_parameters={}, naming=String, build=False, run=True,
                   verbose=True, delay=False):
    '''
    EXAMPLES:

      >>> run_simulation('cat {sim_id}.txt', '{sim_id}.txt', template_text='Hello ${x*10}\n', sweep_parameters={'x': (1, 2, 3)}, verbose=False)
      Hello 10
      Hello 20
      Hello 30
    '''
    if (((template_path is None) and (template_text is None))
        or (not (template_path is None) and not (template_text is None))):
        raise ValueError('Exactly one of `template_path` or `template_text` must be provided.')

    params = single_parameters.copy()

    if build:
        pass

    sim_ids = {}
    keys = list(sweep_parameters.keys())
    values = sweep_parameters.values()

    if sweep_parameters:
        product = itertools.product(*values)
    else:
        product = [params]

    if template_path:
        config = Template(filename=template_path, input_encoding='utf-8')
    else:
        config = Template(text=template_text, input_encoding='utf-8')

    name_gen = naming(length=reduce(operator.mul, [len(value) for value in values], 1))

    processes = []
    for param_set in product:
        sweep_params = params.copy()
        if keys:
            for index, value in enumerate(param_set):
                sweep_params[keys[index]] = value

        sim_id = name_gen.next(keys, sweep_params.values())
        sim_ids[sim_id] = sweep_params

        config_rendered = config.render_unicode(**sweep_params).encode('utf-8', 'replace')
        config_file = open(config_path.format(sim_id=sim_id), 'wb')
        config_file.write(config_rendered)
        config_file.close()
        if run:
            if verbose:
                print("Running simulation {sim_id} with parameters:".format(sim_id=sim_id))
                print('\n'.join('{key}: {param}'.format(key=key, param=param) for key, param in zip(keys, param_set)))
            proc = subprocess.Popen(command.format(sim_id=sim_id), shell=True)
            processes.append(proc)
            if delay:
                time.sleep(delay)

    # Wait until all processes are finished
    for process in processes:
        process.wait()

    return sim_ids
