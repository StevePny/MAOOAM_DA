from simulation_runner import run_simulation, Sequential

import numpy

sim_ids = run_simulation('./maooam_lyap_uncoupled {sim_id}', 'params/params_{sim_id}.nml',
                         template_path='params_template.nml',
                         sweep_parameters={'coupling_motion': numpy.linspace(0.001, 1, 10),
                                           'coupling_thermo': numpy.linspace(0.001, 1, 10)},
                         naming=Sequential(zfill=3), delay=1, run=True)
