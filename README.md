# PyClaw_hbox2D
This is a collection of PyClaw files adapted specifically for using 2D hbox methods to solve 2D shallow water equations (see repo "2D_hbox").

The files included in this repo are only the ones needed to add to and replace the same-named files in the original PyClaw repo, in order to make the examples in "2D_hbox" work:

For diagonal example (located under "diagonal_barrier"):
`pyclaw/src/pyclaw/classic/step2ds.f90`
`pyclaw/src/pyclaw/classic/redist_module.f90`
`pyclaw/src/pyclaw/classic/setup.py`
`pyclaw/src/pyclaw/classic/solver.py`
`pyclaw/src/pyclaw/solver.py`
`pyclaw/src/pyclaw/controller.py`
`pyclaw/src/pyclaw/solution.py`
`pyclaw/src/pyclaw/state.py`

For horizontal barrier example (located under "horizontal_barrier"):
`pyclaw/src/pyclaw/classic/hbox2D_mod.f90`
`pyclaw/src/pyclaw/classic/redist_module.f90`
`pyclaw/src/pyclaw/classic/flux2.f90` 
`pyclaw/src/pyclaw/classic/setup.py` 
`pyclaw/src/pyclaw/classic/step2ds.f90`
`pyclaw/src/pyclaw/classic/solver.py`

Download original PyClaw; go to these directories and either add (if not there) or replace (if there) with these files, one at a time for diagonal and horizontal barrier separately. 

To download original PyClaw, download Clawpack following the instructions in clawpack.org (http://www.clawpack.org/developers.html?highlight=developers%20install#installation-instructions-for-developers). You'll need to download Anaconda3. Once you download Clawpack, you compile by going to `clawpack/pyclaw/src/pyclaw/classic` and run `python setup.py install` after modifying the codes with codes in this repo as discussed above. Then you go to the `build/lib` directory inside the same `classic` directory and copy the `.so` files and paste them to `anaconda3/lib/python3/site-packages/clawpack/pyclaw/classic`. Also, for the diagonal example, you should modify the PyClaw object files in the directory `anaconda3/lib/python3/site-packages/clawpack/pyclaw`.
