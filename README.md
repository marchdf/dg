# Development process

1. Create a branch for the new feature (locally):
	```{bash}
	git checkout -b feature-branch
	```

2. Develop the feature, merging changes often from the develop branch into your feature branch:
	```{bash}
	git commit -m "Developed feature"
	git checkout develop
	git pull                     # fix any identified conflicts between local and remote branches of "develop"
	git checkout feature-branch
	git merge develop            # fix any identified conflicts between "develop" and "feature-branch"
	```

3. Push feature branch to dg1d repository:
	```{bash}
	git push -u origin feature-branch
	```

4. Create a pull request on GitHub. Make sure you ask to merge with develop
