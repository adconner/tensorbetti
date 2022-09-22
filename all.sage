load_attach_mode(load_debug=True)

from itertools import *
from IPython import embed

from sage.libs.singular.option import opt, opt_ctx
from sage.libs.singular.function_factory import ff

attach('ParameterizedVariety.sage')
attach('minors_mod_ideal.sage')
attach('matroid.sage')
attach('tensorbetti.sage')
attach('ideal_util.sage')

# vim: ft=python
