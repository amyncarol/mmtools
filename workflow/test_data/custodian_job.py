from custodian.custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler, FrozenJobErrorHandler, \
        UnconvergedErrorHandler, MeshSymmetryErrorHandler, MaxForceErrorHandler, \
        PotimErrorHandler, NonConvergingErrorHandler, WalltimeHandler 
from custodian.vasp.jobs import VaspJob

vasp_cmd = ['ibrun', '/home1/05018/tg843171/vasp.5.4.4_vtst/bin/vasp_std']
handlers = [FrozenJobErrorHandler(timeout=60)] 
jobs = [VaspJob(vasp_cmd, final=True, suffix="", auto_npar=False)]
c = Custodian(handlers, jobs, max_errors=2)
c.run()

