# Augustus optimiser


This is a script to optimise the extrinsic hints bonus/malus for running augustus.
Eventually it might also offer meta/hyper-parameter optimising for augustus training.


## How not to use this software.

There is a script in the [JAMg](https://github.com/genomecuration/JAMg) pipeline to optimise augustus hints config files.
Unfortunately it only optimises one hint type at a time, so we miss out on optimising combinations of hints.

Augustus already has a pipeline that optimises augustus metaparameters called `optimize_augustus.pl`.
This pipeline seems to work pretty well, but the heuristic that it uses to deal with combinations limits it's parallelism capacity.


## Why should you use this?

Both of the other options use a grid-search strategy, which could miss local patterns in the data if you aren't careful about the range you specify it to search over.
Grid searches are also at a disadvantage when it comes to combinations of hyper-parameters, because a large number of combinations have to be tested to get complete coverage.

Here I'm using a random-search strategy because:

1. It's very simple.
2. It is embarrasingly parallel, so it can be run on a supercomputer quickly with MPI.
3. In the machine learning world it has been shown to find a better solution than grid-searches with fewer iterations (usually).


Another good option might be the bayesian/tree or genetic-algorithm approaches.
I like the look of [DEAP](https://github.com/DEAP/deap).
But let's [KISS](https://en.wikipedia.org/wiki/KISS_principle) for now.


Let's see how it goes!
