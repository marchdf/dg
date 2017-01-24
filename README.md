# Examples

Examples for running the code can be found in `code/examples`.

# Documentation

Documentation is through `doxygen`. Run `make doc` in the `code`
directory to generate the documentation.

# Development process

1. Create a branch for the new feature (locally):
	```{bash}
	git checkout -b feature-branch
	```

2. Develop the feature, merging changes often from the `master` branch into your feature branch:
	```{bash}
	git commit -m "Developed feature"
	git checkout master
	git pull                     # fix any identified conflicts between local and remote branches of "master"
	git checkout feature-branch
	git merge master            # fix any identified conflicts between "master" and "feature-branch"
	```

3. Push feature branch to `dg` repository:
	```{bash}
	git push -u origin feature-branch
	```

4. Create a pull request on GitHub. Make sure you ask to merge with master

# Prerequisites
These are things that need to be installed to compile/run the code:
- [Gmsh](http://www.geuz.org/gmsh/): a mesh generating program
- [GNU Scientific Library](https://www.gnu.org/software/gsl/). Make
  sure you install the dev files (e.g. for Ubuntu: `sudo apt-get
  install libgsl0-dev`
- [BLAS](http://www.netlib.org/blas/)
  and [LAPACK](http://www.netlib.org/lapack/). Make sure you install
  the dev files (e.g. for Ubuntu: `sudo apt-get install libblas-dev
  liblapack-dev`)
- add the following paths to your .bashrc so you can use some nice scripts
  ```
  PATH=$PATH:$HOME/PATHTO/dg/code/scripts
  PYTHONPATH=$PYTHONPATH:$HOME/PATHTO/dg/code/scripts
  ```

# Licensing

See [LICENSE.md](LICENSE.md)

# Author contact

Marc T. Henry de Frahan (marchdf@umich.edu)
