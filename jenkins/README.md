# opm-upscaling jenkins build scripts:

**build-opm-upscaling.sh**:
This is a helper script which contains functions for building,
testing and cloning opm-upscaling and its dependencies.

**build.sh**:
This script will build dependencies, then build opm-upscaling and execute its tests.
It is intended for post-merge builds of the master branch.

**build-pr.sh**:
This script will build dependencies, then build opm-upscaling and execute its tests.
It inspects the $ghbPrBuildComment environmental variable to obtain a pull request
to use for ert, opm-common, opm-parser, opm-material, opm-core, opm-output
and dune-cornerpoint (defaults to master) and then builds $sha1 of opm-upscaling.

It is intended for pre-merge builds of pull requests.

You can optionally specify a given pull request to use for ert, opm-common,
opm-parser, opm-material, opm-core and dune-cornerpoint through the trigger.
The trigger line needs to contain ert=&lt;pull request number&gt; and/or
opm-common=&lt;pull request number&gt; and/or opm-parser=&lt;pull request number&gt;
and/or opm-material=&lt;pull request number&gt; and/or opm-core=&lt;pull request number&gt;
and/or opm-output=&lt;pull request number&gt;
and/or dune-cornerpoint=&lt;pull request number&gt;. asdfasdf
