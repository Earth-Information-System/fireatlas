# How to Release
In this context "releasing" means the folloing things:

* tagging the algorithm with a certain [semantic version](https://semver.org/) (semver for short) 


* building an image off that tag that will be used in some async task runner (currently only DPS) to run the regional algorithm jobs asynchronously.

Most of this can be automated but since [semver](https://semver.org/) is often about
considering if the newest set of changes we are packaging up under a version is backward
compatible it does require a human to choose the version incrementing. 

That said, the `fireatlas` code
isn't a library that others will be using in their code and so that also relieves us of considering
backward incompatible changes. Therefore, we will probably never increment the major release number and only minor or patch
in `<major>.<minor>.<patch>`. Another way of saying that is we `fireatlas` will always be within
`>=1.0.0,<2.0.0` this range

---
## Choose a Version Number

The releaser should look at the [current release tags and versions](https://github.com/Earth-Information-System/fireatlas/tags)
and decide if the minor or patch version should be incremented:

* are all the merged changes in this release just bug fixes? then bump the patch (`1.<minor>.<patch>`) version by one
* did any of the merged changes going out include new features? then bump the minor (`1.<minor>.<patch>`) version by one

## Create a new PR for DPS Jobs

Once the releasaer has a version number, then will need to create a PR that modifies the
algorithm config `algorithm_version` in `./maap_runtime/coordinator/algorithm_config.yaml`:

```yaml
algorithm_description: "coordinator for all regional jobs, preprocess and FireForward steps"
algorithm_version: <NEW VERSION NUMBER HERE>
environment: ubuntu
```

## Merge PR and Manually Release

The release can merge the above PR and then kick off a new release by doing the following:

0. Go to [https://github.com/Earth-Information-System/fireatlas/releases](https://github.com/Earth-Information-System/fireatlas/releases)

1. click "Draft New Release"

2. create a new tag for this release that matches the one chose above

3. click the "Generate release notes"

4. review the release notes and clean up

5. click the "Publish release"


## Release Publish Workflow

The manual step in the last section will kick off an async GH actions workflow that does the following

* uses our version information and builds a python package using `twine`


* publishes the package to PyPI with that version number