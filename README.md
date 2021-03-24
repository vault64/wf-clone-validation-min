# wf-plasmid-assembly

This repository contains a [nextflow](https://www.nextflow.io/) workflow
that can be used to de novo assemble plasmid sequences from Oxford Nanopore data.


## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker of conda is installed.

> See the sections below for installation of these prerequisites in various scenarios.
> It is not required to clone or download the git repository in order to run the workflow.

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-plasmid-assembly --help
```

to see the options for the workflow.

**Workflow outputs**

The primary outputs of the workflow include:

* The `fasta` files for the final polished and intermediate assembly stages
* The `fastq` files for the intermediate datasets constructed by the workflow
* A directory, per sample, of annotations from prokka including a .gff
* an HTML report document detailing the primary findings of the workflow.

### Supported installations and GridION devices

Installation of the software on a GridION can be performed using the command

`sudo apt install ont-nextflow`

This will install a java runtime, Nextflow and docker. If *docker* has not already been
configured the command below can be used to provide user access to the *docker*
services. Please logout of your computer after this command has been typed.

`sudo usermod -aG docker $USER`

### Installation on Ubuntu devices

For hardware running Ubuntu the following instructions should suffice to install
Nextflow and Docker in order to run the workflow.

1. Install a Jva runtime environment (JRE):

   ```sudo apt install default-jre```

2. Download and install Nextflow may be downloaded from https://www.nextflow.io:

   ```curl -s https://get.nextflow.io | bash```

   This will place a `nextflow` binary in the current working directory, you 
   may wish to move this to a location where it is always accessible, e.g:

   ```sudo mv nextflow /usr/local/bin```

3. Install docker and add the current user to the docker group to enable access:

   ```
   sudo apt install docker.io
   sudo usermod -aG docker $USER
   ```

## Running the workflow

The `wf-plasmid-assembly` workflow can be controlled by the following parameters. The `fastq` parameter
is the most important parameter: it is required to identify the location of the
sequence files to be analysed.

**Parameters:**

- `fastq` specifies a path to either a *directory* containing FASTQ files, or a *directory* 
    of such directories. If the former, any fastq files in that single location will be 
    processed together, and treated as a single sample. If the latter, i.e. there are 
    sub-directories, each sub-directory will be expected to contain `fastq` files, and 
    will each be treated as a separate sample. (required)
- `host_reference` FASTA file, reads which map to it are discarded, this is useful for 
   filtering out genomic background (optional)
- `regions_bedfile` BED file used to mask regions within host_reference from filtering (optional)
- `approx_size` Approximate size of the plasmid in base pairs, this need only be a rough estimate (default: 10000)
- `assm_coverage` Try to use this many fold coverage per assembly, note that the pipeline 
   conducts 3 assemblies per sample and reconciles their differences (default: 150)
- `flye_overlap` Sets the min overlap that flye requires of reads (default: 2000)
- `out_dir` the path for the output (default: output)

To run the workflow using Docker containers supply the `-profile standard`
argument to `nextflow run`:

> The command below uses test data available from the [github repository](https://github.com/epi2me-labs/wf-plasmid-assembly/tree/master/test_data)
> It can be obtained with `git clone https://github.com/epi2me-labs/wf-plasmid-assembly`.

```
# run the pipeline with the test data
OUTPUT=output
nextflow run epi2me-labs/wf-plasmid-assembly \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --fastq test_data \
    --out_dir ${OUTPUT}
```

The output of the pipeline will be found in `./output` for the above
example. This directory contains the nextflow working directories alongside
the primary outputs of the pipeline: the intermediate and final assemblies,
which are named `<your_directory>.final.fasta`.

### Running the workflow with Conda

To run the workflow using conda rather than docker, simply replace 

    -profile standard 

with

    -profile conda

in the command above.

### Configuration and tuning

> This section provides some minimal guidance for changing common options, see
> the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for further details.

The default settings for the workflow are described in the configuration file `nextflow.config`
found within the git repository. The default configuration defines an *executor* that will 
use a specified maximum CPU cores (four at the time of writing) and RAM (eight gigabytes).

If the workflow is being run on a device other than a GridION, the available memory and
number of CPUs may be adjusted to the available number of CPU cores. This can be done by
creating a file `my_config.cfg` in the working directory with the following contents:

```
executor {
    $local {
        cpus = 4
        memory = "8 GB"
    }
}
```

and running the workflow providing the `-c` (config) option, e.g.:

```
# run the pipeline with custom configuration
nextflow run epi2me-labs/wf-plasmid-assembly \
    -c my_config.cfg \
    ...
```

The contents of the `my_config.cfg` file will override the contents of the default
configuration file. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html)
for more information concerning customized configuration.

**Using a fixed conda environment**

By default, Nextflow will attempt to create a fresh conda environment for any new
analysis (for reasons of reproducibility). This may be undesirable if many analyses
are being run. To avoid the situation a fixed conda environment can be used for all
analyses by creating a custom config with the following stanza:

```
profiles {
    // profile using conda environments rather than docker
    // containers
    fixed_conda {
        docker {
            enabled = false
        }
        process {
            withLabel:artic {
                conda = "/path/to/my/conda/environment"
            }
            shell = ['/bin/bash', '-euo', 'pipefail']
        }
    }
}
```

and running nextflow by setting the profile to `fixed_conda`:

```
nextflow run epi2me-labs/wf-plasmid-assembly \
    -c my_config.cfg \
    -profile fixed_conda \
    ...
```

## Updating the workflow

Periodically when running the workflow, users may find that a message is displayed
indicating that an update to the workflow is available.

To update the workflow simply run:

    nextflow pull epi2me-labs/wf-plasmid-assembly

## Building the docker container from source

The docker image used for running the `wf-plasmid-assembly` workflow is available on
[dockerhub](https://hub.docker.com/repository/docker/ontresearch/wf-plasmid-assembly).
The image is built from the Dockerfile present in the git repository. Users
wishing to modify and build the image can do so with:

```
CONTAINER_TAG=ontresearch/wf-plasmid-assembly:latest

git clone https://github.com/epi2me-labs/wf-plasmid-assembly
cd wf-plasmid-assembly

docker build \
    -t ${CONTAINER_TAG} -f Dockerfile \
    --build-arg BASEIMAGE=ontresearch/base-workflow-image:v0.1.0 \
    .
```

In order to run the workflow with this new image it is required to give
`nextflow` the `--wfversion` parameter:

```
nextflow run epi2me-labs/wf-plasmid-assembly \
    --wfversion latest
```

## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)