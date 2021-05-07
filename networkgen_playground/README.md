# Network generator playground #
This folder contains development for the implementation of the interactive
network to be used in qfepweb. Current implementation uses [the community version of visjs](https://github.com/visjs-community/visjs-network) which allows for high-level network generation with scalability to >1000 nodes. 

This environment has been written in a way to facilitate development; nearing the point of django integration data flow design likely has to be adjusted. Currently the workflow works by starting from `./input/`, using `./graphgen.py` (which imports `./helper_functions.py`) to write data to `./data/`. `./networkgen.html` ports data from `./data/` and broadcasts an interactive network while writing output to `./output/`.

Also one can follow issue #1 to follow the discussion and contribute.  

https://github.com/GPCR-ModSim/qfepweb/issues/1  

Developers:

- Jenke Scheen

  

**Previous version:**

The notebook can be perused using mybinder.org (Just click on the binder icon at the bottom).  



[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/GPCR-ModSim/qfepweb/networkgen?filepath=networkgen_playground%2Fgraph_gen_tests.ipynb)
