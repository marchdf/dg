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


# Licensing

See [LICENSE.md](LICENSE.md)

# Author contact

Marc T. Henry de Frahan (marchdf@umich.edu)
