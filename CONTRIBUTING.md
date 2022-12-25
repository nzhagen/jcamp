# Contributing to `jcamp`


## Submit a patch / proposal

Your contributions and hints are welcome!

Please [open an issue](https://github.com/nzhagen/jcamp/issues) if you want to discuss problems
or discuss improvements.

You code contributions are welcome, too.
*Pull requests* via github are the preferred approach for these.
Please take a look at the
[github documentation](https://docs.github.com/en/get-started/quickstart/contributing-to-projects)
if you are new to this platform.


## Publish a new release

1. update `__version__` in `jcamp.py`
1. update `version` in `setup.py`
1. create a release commit: `git commit -m "Release v0.1.2" jcamp.py setup.py`
1. attach a tag to the release commit: `git tag v0.1.2`
1. push the commit and its tag: `git push --follow-tags`
1. publish the release to [PyPI](https://pypi.org/) (TODO: describe the procedure)
