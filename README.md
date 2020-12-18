# PyClaw_hbox2D
This is a collection of PyClaw files adapted specifically for using 2D hbox methods to solve 2D shallow water equations (see repo "2D_hbox").

The files included in this repo are only the ones needed to add to and replace the same-named files in the original PyClaw repo, in order to make the examples in "2D_hbox" work:

For diagonal example:
`pyclaw/src/pyclaw/classic/step2ds.f90`
`pyclaw/src/pyclaw/classic/redist_module.f90`
`pyclaw/src/pyclaw/classic/setup.py`
`pyclaw/src/pyclaw/classic/solver.py`
`pyclaw/src/pyclaw/solver.py`
`pyclaw/src/pyclaw/controller.py`
`pyclaw/src/pyclaw/solution.py`
`pyclaw/src/pyclaw/state.py`

For horizontal barrier example:
`pyclaw/src/pyclaw/classic/hbox2D_mod.f90`
`pyclaw/src/pyclaw/classic/flux2.f90` 
`pyclaw/src/pyclaw/classic/setup.py` 
`pyclaw/src/pyclaw/classic/step2ds.f90`
`pyclaw/src/pyclaw/classic/solver.py`

Download original PyClaw; go to these directories and either add (if not there) or replace (if there) with these files, one at a time for diagonal and horizontal barrier separately. 
