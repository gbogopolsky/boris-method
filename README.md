# The Boris-Buneman particle pusher

A simple Python implementation of the Boris-Buneman particle pusher.

We then study the trajectory of a particle in ExB configuration and the magnetic configuration of the
Tore-Supra tokamak in Cadarache, France.

This work was realised in collaboration with Lucas Fuster, in the context of a practical session
of the Physics of Plasma and Fusion M. Sc. at Université Paris-Saclay and Sorbonne Université.

## Running the scripts

Use the `conda` package manager and the requirements to create the corresponding Python
environment.
To do so, use the following instructions:
```bash
conda create --name boris --file requirements.txt
source activate boris
python3 the_file.py
```

You can execute three different scripts:
- `main.py` for the ExB configuration
- `error.py` to check the stability of the method vs. the timestep
- `tokamak.py` to compute the trajectory of a particle in Tore Supra.

## Thanks

I would like to thank Andrea Ciardi from the [LERMA](https://lerma.obspm.fr/) lab at OBSPM for putting together this practical, of which
you may find the files [here](https://sites.google.com/site/andreaciardihomepage/teaching?authuser=0).
