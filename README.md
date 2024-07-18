## **Standardizing and harmonizing NGS analysis workflows**

Materials for NGS Harmonization Workshop at GCB 2024

It is planned to be 3 hours of an overview of standardization and harmonizing NGS analysis strategies in GHGA. We will explore how FAIR principles enable the standardization and harmonization of nf-core-based NGS analysis workflows within GHGA. We will  demonstrate the adaptability of nf-core workflows and discuss the importance of standardization of workflows. Finally, we will show how to make workflows scalable, robust, and automated using a small subset of a public dataset. 

### Preliminary Schedule

|Time|Topic|
|:---|:---|
|9:00am  - 9:10am| Introduction to the tutorial: What is GHGA? What are our workflow objectives? What is FAIR data|
|9:10am  - 9:30am| Reproducibility, adaptability, and portability of Workflows|
|9:30am  - 10:15am| Hands-on part 1: Docker, Github, introduction to nextflow and nf-core tools |
|10:15am - 10:30am| Break|
|10:30am - 11:30am| Hands-on part 2: Setting up a simple nextflow workflow using nf-core tools|
|11:30am - 11:50am| How to keep reproducibility 
|11:50am - 12:00am| Summary 


### Learning Objectives for Tutorial

- FAIR principles and their relevance for workflow standardization and harmonization
- The adaptability, portability, and scalability of workflows within nf-core 
- Automating robust workflows to ensure reproducibility and the highest quality of code

### Overview

We will walk through the following steps:
1. Learning how to collaborate using GitHub.
2. Constructing an efficient workstation through Visual Studio Code. 
3. Build and use Docker containers to encapsulate software dependencies. 
4. Set up a development environment to run Nextflow
5. Overview of nf-core tools, pipelines, and modules.
6. Contraction of a simple variant calling pipeline using bwa-mem and bcftools using nf-core modules and templates.
7. Dealing with config files and running the pipeline. 
8. Exploring the results
9. Adapting a Continuous Integration environment using github actions. 

### Requirements
Please have the following software and user accounts ready on the day of the workshop.
- [GitHub Account](https://github.com/)
- [Docker Hub Account](https://hub.docker.com/signup)
- [Visual Studio Code](https://code.visualstudio.com/) or your favorite code editor. 

A preconfigured Nextflow development environment is available using Gitpod. To run **Gitpod**:

- Click the following URL: https://gitpod.io/#https://github.com/nextflow-io/training
  -- This is nextflows GitHub repository URL, prefixed with https://gitpod.io/#
- Log in to your GitHub account (and allow authorization).
- Once you have signed in, Gitpod should load (skip prebuild if asked).
- _If you decided to use the Gitpod environment, you don`t need to install anything into your local computer! But still, you would need to have github and docker accounts_

To follow the workshop on your computer, you will need the following software and files:
- Bash
- [Java11](https://www.oracle.com/java/technologies/downloads/) (or later, up to 18)
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- [nf-core](https://nf-co.re/) - [nf-core tools on Github](https://github.com/nf-core/tools)
- [Docker](https://www.oracle.com/java/technologies/downloads/)
- [GitHub CLI](https://cli.github.com/)
- [Data](https://github.com/GHGA-Training/gcb_ngs_harmonisation_workshop/tree/main/data) 

## Part 1 

### 1. What is **GitHub**, and how we can use it?

The most common way of providing stable versions for software is through git systems. Github, gitlab or bitbucket are all valuable services in software development. They provide code hosting/storage platform for version control and collaboration. It lets you and others work on the same project from anywhere. Using git systems one can create different branches (or copies) of the stable working version dependent or independent of one another. 

GitHub is one of the most commonly used platform for collaborative software development and version control using git. We will be using the following repositories through this workshop. If you want, you can fork [GHGA-Training/gcb_ngs_harmonisation_workshop](https://github.com/GHGA-Training/gcb_ngs_harmonisation_workshop) or follow it from the original documentation.

Go to the websites and fork the following repositories:
- https://github.com/GHGA-Training/gcb_ngs_harmonisation_workshop

or

```bash
gh auth login
gh repo fork https://github.com/GHGA-Training/gcb_ngs_harmonisation_workshop.git --clone
```

Whenever you fork the page, you can take any actions on that as you wish! 

Now, lets clone the forked gcb_ngs_harmonisation_workshop and start working on it.

```bash
git clone https://github.com/<yourgithubname>/gcb_ngs_harmonisation_workshop
git checkout -b dev
```

Usually, working on branches is advised to keep master environment as safe as possible. So, lets create _testbranch_ branch:

```bash
git checkout -b testbranch
```

- You can switch back to master:

```bash
git checkout master
```

- You can delete the branches like: 

```bash
git checkout -d testbranch
```

Don't forget that, unless you push the branch to your remote repository, the changes are only available to you.

- You can push the branch to your remote like:

```bash
git push origin testbranch
```

Exercise:

- Now, checkout to a new branch named as '_dev_' and start working on it.  

> _*Note*: if you don't want to fork the repo, pull it in your local computer work on it and push it into your own github repo!_ 

### 2.  Creating your own VS-Code workstation:

There are plenty of good editors but here we will be using VS-Code. If you have any other preference just go with it! The idea is to be able to connect multi-services in an environment creating a functional development strategy for lining, debugging,  editing, and pushing your code. You can add, VS-code extensions like Nextflow, nf-core tools, github, python, R and many more. 

Simply, you can open _gcb_ngs_harmonisation_workshop_ using VS-code and start working on it!

- Open VS-code

- if you did not clone the repository:
    - Press "Clone Git Repository"
- Else:
    - Open the directory where you clone the repository

- Checkout to '_dev_' branch 
    - Press '_master_' bottom left of the screen
    - Write dev and press enter.
  
### 3.  How to build and use **Docker** containers:

- Containers are configurable virtualization technology that allows the packaging and distribution of pipelines  in a self-contained and platform-independent manner. By this way, installed software using package managers can be packed into an image with its corresponding dependencies. An image can be considered as  a file (or set of files) that contains the application with the code and all its dependencies, libraries, etc. You can copy images around, upload them into registries, download them, and re-use them. Commonly used software for containerization in bioinformatics includes Docker, Singularity and podman. 

- Using software containers is crucial to make our pipelines portable. In this workshop, we will be using docker and singularity in order to dive in better, we will be constructing our own bcftools container as an example.  
  
**Docker**

```bash
docker login
```

- The pull command lets you download a Docker image without running it. For example

```bash
docker pull ubuntu:focal
```

 - Launching a BASH shell in the container allows you to operate in an interactive mode in the containerized operating system. For example

```bash
docker run -ti ubuntu:focal
```
- Now, you will be inside the ubuntu:focal image. But, in your local system, another layer with root@<containerid> should have been generated. 
  
- Lets build an image to download and install *bcftools*

```
# Update the package inside the container
apt-get update
# Install the tools we need to download and compile Samtools
apt-get install -y bzip2 g++ libbz2-dev libcurl4-openssl-dev liblzma-dev make ncurses-dev wget zlib1g-dev
# Download Bcftools
wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2

# Unpack the archive
tar jxf bcftools-1.17.tar.bz2
cd bcftools-1.17
# Compile the code
make
# Install the resulting binaries
make install
#check if everything is properly installed
bcftools --version
#If complete exit docker
exit
```
root@<containerid> is now exit. And you are back to your local computer. 

- Now, lets commit (save) the <containerid> to our docker registry and upload it there.
(An example image name could be: kubran/bcftools:v1.17)

Prepare your own credentials: 

```bash
docker commit <containerid> <docker_name>/imagename:version_tag
```

- Push the image to Dockerhub registry when the commit is done:

```bash
docker push <docker_name>/imagename:version_tag
```

**Dockerfiles**

- The repeatability of the containerization can be enabled through Dockerfile. A Dockerfile is a text document that contains all the commands a user could call on the command line to assemble an image. It has a big advantage over manuel containerization since the file can be versioned, maintained, edited and shared with others easily. 

- Now, lets use a Dockerfile to create another image: 
    - Save this Dockerfile to your environment (don't use any extension and name the file as Dockerfile):

```
# Use an appropriate base image
FROM ubuntu:20.04

# Set the maintainer label
LABEL maintainer="yourname@example.com"

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    libbz2-dev \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    gcc \
    make \
    && rm -rf /var/lib/apt/lists/*

# Set environment variables
ENV BCFTOOLS_VERSION=1.17

# Download and install BCFtools
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -xvjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    make && \
    make install && \
    cd .. && \
    rm -rf bcftools-${BCFTOOLS_VERSION} bcftools-${BCFTOOLS_VERSION}.tar.bz2

```

- Run the following to build the image from Dockerfile:

```bash
docker build --platform linux/amd64 -t bcftools:1.17 .
```

- Tag the image to Dockerhub registry

```bash
docker tag  bcftools:1.17 <docker_name>/bcftools:1.17
```

- Push the image to Dockerhub registry

```bash
docker push <docker_name>/bcftools:1.17
```

**Singularity**

- Singularity is another containerization platform that you can use instead of docker. Docker images can quickly be converted into singularity. It is also possible to use docker-converted singularity images directly by nextflow as follows:

singularity_image = docker://<docker_name>/imagename:version_tag

### 4.  Getting familiar with **Nextflow** and **nf-core**

The new generation way of pipelines are  Workflow management systems providing all components of standard analysis. They increase transparency, enable long-term sustainability, and aid in achieving findable, accessible, interoperable, and reusable computational analysis and they are FAIR-approved! Snakemake, Galaxy, cornwell, pegasus and Nextflow are all open source and free-to-use examples! Workflow managers use containers to run software providing portability and stability. They provide scalability: They automatically manage the scheduling enabling the effective utilization of available resources. Required reference files can be pulled through cloud registries ensuring saving space locally. They are able to cache intermediate results avoiding recalculations saving significant time and computing resources which is a key advantage of workflow managers. Project repositories are often in git systems allowing collaboration and feedback. When the execution is finished, They provide the option to generate execution reports with detailed information about the run. 

Nextflow is growing fast and has long-term support available from Seqera Labs. It has been developed since 2013 by the same team. Nf-core is a community effort to collect a curated set of analysis pipelines built using Nextflow. The effort itself is not just restricted to pipelines, they also provide sub-workflows and modules (single processes running tools). Their modular design is very standard and also flexible for multi-purposing.  They also provide helper tools to create template pipelines, install modules, and create test cases. 


```bash
nextflow info
```

**nf-core** is a community effort to optimize and collect a set of nextflow pipelines. Currently, they have **108** pipelines!

- let's explore **nf-core** tool:

```bash
nf-core --help
```

- Listing available pipelines:

```bash
nf-core list
```

Running nf-core pipelines are very easy! Most of the resources that nf-core uses are deposit on the cloud, meaning that you might not need to deal with pipeline and reference set-up. Lets perform an a test run using their famous [rnaseq pipeline](https://nf-co.re/rnaseq/3.14.0) (or you can choose to run another pipeline in their list). 

- [_test_ profile](https://github.com/nf-core/rnaseq/blob/3.14.0/conf/test.config) includes necessary input, output and parameter information to perform an example run

```bash
# Launch the RNAseq pipeline
nextflow run nf-core/rnaseq -profile test,docker
```

**nf-core** also provides ready to use atomic processes for pipelines which are called as _modules_. Their modules library is very wide, they have **1273** modules.

- Listing available modules through nf-core:

```bash
nf-core modules list remote
```

We will play more with nf-core modules in the second part. 


## Part 2 

## Development of a simple alignment and variant calling pipeline using nf-core tools


### 1.  Using nf-core modules to build a simple variant calling pipeline

> _*Note*: Creating a new pipeline using nf-core commends is well explained [here](https://nf-co.re/docs/tutorials/nextflow_training/creating_with_nf-core)_

- We will use nf-core template to create our first nextflow pipeline. Run the create command:

```bash
nf-core create
```

I used "mydemo" to name my pipeline. nf-core automatically adds nf-core prefix if you accept default settings. You can follow the instructions and interactive prompts.

- Lets move into the new directory created with template:

```bash
cd nf-core-mydemo
```

- nf-core template comes with fully ready pipeline, before we go further be the command also initiated a git repository for you:

```bash
git status
git branch
```

It is a common practice to work on _dev_ or _master_ branches, you should not touch _TEMPLATE_. nf-core has a automated template synchronization and it works through with template branch.

- Since the template created a full ready pipeline, we can actually run it! It has *FASTQC* and *MultiQC* modules installed with a small test config which contains small sets of reads.

```bash
cd ../
nextflow run nf-core-mydemo/ -profile test,docker --outdir test_results
```

- Check test_results/ output directory to study fastqc outputs. 

---



Lets start building our own pipeline now:

- We will perform _bwa-mem_ alignment, which requires indexed fasta genome from _bwa-index_. Thus we need bwa-mem and bwa-index modules. Luckily, nf-core modules library provides bwa modules readily available. Lets install bwa/mem and bwa/index modules:

```bash
cd nf-core-mydemo
nf-core modules install bwa/mem
nf-core modules install bwa/index
```

Now both modules should be located in **modules/nf-core/bwa** directory. Check main.nf files to understand input and output structures.

BWA_INDEX module needs a tuple fasta file and outputs a tuple bwa/ directory containing fasta indexes. Those two files are both attached to the same meta. Which will enable us to track files easily. 

```nextflow
process BWA_INDEX {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.18--he4a0461_0' :
        'biocontainers/bwa:0.7.18--he4a0461_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(bwa) , emit: index
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    def args   = task.ext.args ?: ''
    """
    mkdir bwa
    bwa \\
        index \\
        $args \\
        -p bwa/${prefix} \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    mkdir bwa

    touch bwa/${prefix}.amb
    touch bwa/${prefix}.ann
    touch bwa/${prefix}.bwt
    touch bwa/${prefix}.pac
    touch bwa/${prefix}.sa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}

```

- In order to **include** those modules to our pipeline, we will add descriptions/paths to the workflow. 

Add two lines of code to workflows/mydemo.nf

```nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

...
include { BWA_MEM                } from '../modules/nf-core/bwa/mem/main'
include { BWA_INDEX              } from '../modules/nf-core/bwa/index/main'
...
```

- To perform alignment using _bwa-mem,_ we need the reference fasta file to be indexed. 

>_*Note*: The template includes igenome.config template ready to use. Therefore, we can directly use one of the provided fasta files readily. We will use _bwa_index_ tool to index fasta file. [IGenomes](https://emea.support.illumina.com/sequencing/sequencing_software/igenome.html) is a source providing collections of references and annotations supported by AWS. igenome.config includes parameters for the available sources for the pipeline._ 

The template comes ready to use fasta channel. Please take a look into main.nf file. There _getGenomeAttribute_ function enables initiation of the channel from parameters.

- Now, cut the following part in main.nf and put it below _nextflow.enable.dsl = 2_ line. 

!! This part is only necessary due to a current bug in the template!!

```nextflow
nextflow.enable.dsl = 2
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_variantbenchmarking_pipeline'
params.fasta = getGenomeAttribute('fasta')

```

- Lets create _fasta_ now.  Open workflows/mydemo.nf and paste the following lines:

```nextflow

workflow MYDEMO {

    take: 
        samplesheet // channel: samplesheet read in from --input

    main:
...
    // check mandatory parameters
    println(params.fasta)
    ch_fasta       = Channel.fromPath(params.fasta, checkIfExists: true).map{ it -> tuple([id: it[0].getSimpleName()], it) }.collect()

...

}

Now, 

- and lets have a look on ch_fasta, type this to mydemo.nf and make a test run.

```nextflow

ch_fasta.view()
```

Run the pipeline:

```bash
cd ../
nextflow run nf-core-mydemo/ -profile test,docker --outdir test_results --fasta "mydata/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
```

> _*Note*: What is really cool and helpful to use the **-resume** tag to resume previously finished jobs! Moreover, don't forget to check out .nextflow.log files in case of an error. All of the runs will be saved into _work_ directory._


- Now, we are ready to add our first module! Checking out BWA_INDEX, the only input file is a fasta file to be able to create bwa index directory. Before you go through the process take a few minutes to read about modules [here](https://www.nextflow.io/docs/latest/module.html). 

The output index directory will be saved into ch_index channel for our further usage. 

```nextflow
    //
    // MODULE: BWA_INDEX
    //
    BWA_INDEX(
        ch_fasta
    )
    ch_index = BWA_INDEX.out.index
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
```

Run the pipeline to see the results:
```
nextflow run nf-core-mydemo/ -profile test,docker --outdir test_results --fasta "mydata/ggal_1_48850000_49020000.Ggal71.500bpflank.fa" -resume
```

> **Exercise:** Check BWA_INDEX output directory.   

Next task is to add BWA_MEM to align our reads, but first lets edit samplesheet.csv so that it will include our own samples:

- This pipeline comes with nf-core sub-workflow in order to check the input files and create an input channel with them (subworkflows/nf-core/utils_nfvalidation_plugin/main.nf). It automatically checks and validates the header, sample names, and sample directories. nf-validation plug-in uses assets/schema_input.json to automatically detect csv header from samplesheet.csv and makes it easy to use as input_ch. assets/schema_input.json file in this project tuned to use sample, fastq_1 and fastq2 for example, since we will be using the same format in our analysis we are not going to change schema_json. Please read more [here](https://nextflow-io.github.io/nf-schema/latest/) about nf-schema structure. 

- Lets create our own samplesheet.csv file using the test samples given:
     - create mysamplesheet.csv file and save it under assets/ directory: 


```console
sample,fastq_1,fastq_2
gut,mydata/ggal_gut_1.fq.gz,data/ggal_gut_2.fq.gz
liver,mydata/ggal_liver_1.fq.gz,data/ggal_liver_2.fq.gz
```

!! Make sure to give the correct path for your files (instructions were given already for sample files)


When we do test runs, we were using test.config. That file is readily available for template pipeline describing the parameters needed to run this pipeline involving _--input_ and _--genome parameters.

- Now, lets create our own config file:
    - Open a text file and create a config file as follows and name to mytest.config


``` nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow config file for running minimal tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Defines input files and everything required to run a fast and simple pipeline test.

    Use as follows:
        nextflow run nf-core/nf-core-mydemo -profile mytest,<docker/singularity> --outdir <OUTDIR>

----------------------------------------------------------------------------------------
*/

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Limit resources so that this can run on GitHub Actions
    max_cpus   = 2
    max_memory = '6.GB'
    max_time   = '6.h'

    // Input data
    input  = 'assets/mysamplesheet.csv' 

    // Genome references
    //genome = 'hg38'
    fasta = "mydata/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
}
```

Now we will be adding BWA_MEM alignment module:

BWA_MEM module requires fastq reads which were already provided with ch_samplesheet. We prepared ch_index channel in the previous step using _BWA_INDEX_ module.

"_sort_" statement is a value in order to activate samtools sort for BAM file.

- add BWA_MEM module to mydemo.nf : 

```nextflow
    //
    // MODULE: BWA_MEM
    //
    BWA_MEM(
        ch_samplesheet,
        ch_index,
        fasta,
        "true"
    )
    ch_bam = BWA_MEM.out.bam
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)
```

Run the pipeline to see the results:

```bash
nextflow run nf-core-mydemo/ -profile docker -c conf/mytest.config --outdir test_results -resume
```

> Warning: Be awaire that we are no longer using the test profile and switched into our own config. To be able to run it through profile, we should add it as profile in nextflow.config:

```nextflow
    test      { includeConfig 'conf/test.config'      }
    test_full { includeConfig 'conf/test_full.config' }
    mytest    { includeConfig 'conf/mytest.config'    }
```

Final task will be calling variants using the bam files generated with BWA_MEM. We will use bcftools mpileup together with bcftools call and bcftools view. Let's create our first module using bcftools container that we created!

> Tip: If you think this task is too much for you, you can also use readly available BCFTOOLS_MPILEUP module from nf-core! Just follow the same steps that we did for BWA_INDEX and BWA_MEM!_


- A draft BCFTOOLS_MPILEUP module will look like this: We will need to define input and output files together with bcftools commands that will process variant calling. 

```nextflow
process BCFTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'singularity_container':
        'docker_container' }"

    input:
    INPUT_FILES

    output:
    OUTPUT_FILES
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    CMD

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CMD_VERSION
    END_VERSIONS
    """
}

```
Open bcftools_mpileup.nf and place under modules/local folder, and copy the above draft BCFTOOLS_MPILEUP module into open bcftools_mpileup.nf. We will need to include _INPUT_FILES_, _OUTPUT_FILES_, _CMD_ and _CMD_VERSION_ sections in order to make this module compleate as follows: 

Now, we need to construct the simple CMD for variant calling: 

```nextflow
    bcftools \\
        mpileup \\
        --fasta-ref $fasta \\
        -Oz \\
        $bam \\
        | bcftools call --output-type v -mv -Oz \\
        | bcftools view --output-file ${prefix}.vcf.gz --output-type z

    tabix -p vcf -f ${prefix}.vcf.gz
    bcftools stats ${prefix}.vcf.gz > ${prefix}.bcftools_stats.txt

```
_bctools mpileup_ will produce a mpileup file including genotypes, then we will use _bcftools call_ to actually call variant sites and save as gzipped vcf file. As a plus, we will use _bcftools stats_ to examine the number of variants. 

We will need the alignment bam file and reference fasta file to run _bcftools mpileup_ and this argument will produce an indexed vcf.gz file and a statistic file in txt format. Therefore, input and output definitions will be: 

```nextflow
    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path (fasta)

    output:
    tuple val(meta), path("*vcf.gz")     , emit: vcf
    tuple val(meta), path("*vcf.gz.tbi") , emit: tbi
    tuple val(meta), path("*stats.txt")  , emit: stats
    path  "versions.yml"                 , emit: versions

```

 We will also emit versions file to keep track of bcftool versioning.


```nextflow
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
```

 We can add either our docker container (created in the 3rt step) or search for the proper bcftool version on[ biocontainers registry](https://quay.io/repository/biocontainers/bcftools?tab=tags). 

```nextflow
    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'biocontainers/bcftools:1.17--haef29d1_0' }"
```

 The final module should like this:


```nextflow
process BCFTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path (fasta)

    output:
    tuple val(meta), path("*vcf.gz")     , emit: vcf
    tuple val(meta), path("*vcf.gz.tbi") , emit: tbi
    tuple val(meta), path("*stats.txt")  , emit: stats
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools \\
        mpileup \\
        --fasta-ref $fasta \\
        -Oz \\
        $bam \\
        | bcftools call --output-type v -mv -Oz \\
        | bcftools view --output-file ${prefix}.vcf.gz --output-type z

    tabix -p vcf -f ${prefix}.vcf.gz

    bcftools stats ${prefix}.vcf.gz > ${prefix}.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
```

- Since we finished BCFTOOLS_MPILEUP module, we can include it in the workflow (mydemo.nf):

```nextflow
//
// MODULE: Installed directly locally
//
include { BCFTOOLS_MPILEUP            } from '../modules/local/bcftools_mpileup.nf'

```
  
- Let's connect BCFTOOLS_MPILEUP module to the mapped reads output from BWA_MEM module:

```nextflow
    //
    // MODULE: BCFTOOLS_MPILEUP
    //
    BCFTOOLS_MPILEUP(
        ch_bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)
```

One of the other configurations that we need to be careful of is setting a mirror to your container settings. In the tempalate pipeline they set _quay.io_ as a mirror for the registries by default. But if you want to use your own docker container generated in the 3rd step you should delete the default settings. 

```nextflow
apptainer.registry   = 'quay.io'
docker.registry      = 'quay.io'
podman.registry      = 'quay.io'
singularity.registry = 'quay.io'
```
to 

```nextflow
apptainer.registry   = ''
docker.registry      = ''
podman.registry      = ''
singularity.registry = ''
```

If you decide to use ready containers from biocontainers default settings can remain. 


### 2.  Our first full-functioning pipeline is ready! and we can directly run it!

It is actually that easy to create such a functional workflow! 

Here, we will try to run our workflow using the config file we created (mytest), and input file we prepared (mysamplesheet.csv) and docker environment. But you can also run it with singularity. You can also define another output directory else than _test_results_. 

```bash
nextflow run nf-core-mydemo/ -profile mytest,docker --outdir test_results -resume

```

In the first run of this workflow, container images and genome would need to be pulled (all automatic) which is why the initial run will take a bit longer than the following runs!  

What is really cool about nf-core tools is that they could lint your project and tell what it is missing. Lets lint our current pipeline:

```bash
nf-core lint
```

This line of code will print out some warnings about _TODOs_ automatically involved with template pipeline. You can go ahead and try to fullfill the requirements yourself. You can write a nice _README.md_ describing how the pipeline works,or add some extra parameters using _nexftlow.config_ and nextflow_schema.json_. 


### 3. Defining an output structure

If you analyzed the test_results directory already, you would see that it automatically created _bwa_ and _bcftools_ directories for you. Those directories actually are initials of our module names. In order to shape our output structure better lets investigate conf/modules.config which contains arguments per module options:

```nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Config file for defining DSL2 per module options and publishing paths
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Available keys to override module options:
        ext.args   = Additional arguments appended to command in module.
        ext.args2  = Second set of arguments appended to command in module (multi-tool modules).
        ext.args3  = Third set of arguments appended to command in module (multi-tool modules).
        ext.prefix = File name prefix for output files.
----------------------------------------------------------------------------------------
*/

process {

    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    withName: FASTQC {
        ext.args = '--quiet'
    }

    withName: 'MULTIQC' {
        ext.args   = { params.multiqc_title ? "--title \"$params.multiqc_title\"" : '' }
        publishDir = [
            path: { "${params.outdir}/multiqc" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }

}
```

- Lets add shape our outputs better for *BWA_MEM* and *BCFTOOLS_MPILEUP*. Add the following lines to modules.config

```nextflow
...

    withName: 'BWA_MEM' {
        publishDir = [
            path: { "${params.outdir}/${meta.id}/alignment" },
            pattern: "*{.bam}",
            mode: params.publish_dir_mode
        ]
    }
    withName: 'BCFTOOLS_MPILEUP' {
        publishDir = [
            path: { "${params.outdir}/${meta.id}/variant_calling" },
            pattern: "*{.vcf.gz,tbi,stats}",
            mode: params.publish_dir_mode
        ]
    }    
...
```

Please read more from [here](https://www.nextflow.io/docs/latest/config.html) to better learn and understand about nextflow configurations.

### 3. Optimizing resources

nf-core modules uses the average amount of resources to run the process inside. Therefore, they are being optimized into the tool and the input size (average size). Running through some of the modules might require extensive amount of memory and cpu like the case we have in BWA_MEM which uses 'process_high'. In order to run an small test sample, we won't need high processing but can use small amounts. Let's change it. 

We can redefine sources (cpu=4 and memory=6GB) in mytest.config by adding these lines:

```nextflow
process {
   withName:'BWA_MEM'{
        cpus            = { check_max( 2 * task.attempt, 'cpus' ) }
        memory          = { check_max( 6.GB * task.attempt, 'memory' ) }
   }
}
```

> _*Note*: If the resource is low for the running process, task.attempt will automatically be multiplied by the resource and task will be rerun!_

When you create your own pipelines, you can also edit base.config where process tags are being defines.

### 4. Analyzing the results

- _Provenance_ describes the history of a computational experiment and it is a very important aspect of workflows. Minor changes, usually in the workflow environment, software versions, parameter settings, and reference annotations, may or may not change the results of each workflow run but the major challenge is to know and record the exact environment you run your analysis. Data provenance describes this trail of methods, versions, and arguments that were used to generate a set of files. In short, it is the history of a computational environment you used for the analysis.

- nf-core pipelines are armed through many provenance tools with the construction of _work_ directory with _temp_ results of each run, saving workflow execution history through .nextflow.log files as well as _CUSTOM_DUMPSOFTWAREVERSIONS_ module and MultiQC tools. 

- Collection of versions is a vital process in order to keep track of software history. In nf-core pipelines, each tool version is collected in a channel and then processed using _CUSTOM_DUMPSOFTWAREVERSIONS_ module and represented through MultiQC tool.

- Used data, produced inputs and outputs, artifacts, information about the execution like consumed memory and CPU, processed time saved and represented in final reports of nf-core tools. 

- MultiQC tool also aggregates logs and reports from the analysis. In our analysis, FASTQC analysis was already included. In this example file, you can both see the FASTQC report and the software versions together with the workflow summary.


### 5. Committing changes to github repo, and Continuous Integration test environment  

Another advantage of using github is easy debugging. Testing the code and building can be automized. Continuous Integration (CI) tests refer to the build and unit testing stages of the software release process and are curious for automatic software development. Every revision that is committed triggers an automated build and test. When it is reviewed and approved, it can be merged safely. The aim is to protect the main branch from errors and unintended changes. CI tests are mandatory for nf-core pipelines and luckily in our example nf-core-mydemo also comes with an environment that we can actually work on! In order to activate automatic tests:

- Activate actions tab on github (on your github repo nf-core-mydemo)
- under .github/workflows/ci.yml

Change this line

```yml
    # Only run on push if this is the nf-core dev branch (merged PRs)
    if: "${{ github.event_name != 'push' || (github.event_name == 'push' && github.repository == 'nf-core/mydemo') }}"
```
to
```yml
    # Only run on push
    if: "${{ github.event_name != 'push' || (github.event_name == 'push' && github.repository == '<your_github_name>/nf-core-mydemo') }}"
```

We should also change the test profile to ours: 

change this line

```yml
      - name: Run pipeline with test data
        # TODO nf-core: You can customise CI pipeline run tests as required
        # For example: adding multiple test runs with different parameters
        # Remember that you can parallelise this by using strategy.matrix
        run: |
          nextflow run ${GITHUB_WORKSPACE} -profile test,docker --outdir ./results
```

to 

```yml
      - name: Run pipeline with test data
        # TODO nf-core: You can customise CI pipeline run tests as required
        # For example: adding multiple test runs with different parameters
        # Remember that you can parallelise this by using strategy.matrix
        run: |
          nextflow run ${GITHUB_WORKSPACE} -profile mytest,docker --outdir ./results
```

Now whenever we commit changes to _dev_ branch, a small test pipeline will be triggered in order to check if the workflow/code works smoothly.  

    
The last thing we need to do is to commit the changes we have in our repository and push into _dev_ branch that we were working in. 

```bash
git add .
git commit -a -m "commit" (commit message)
```

It is also possible to commit changes using VS-Code if you installed github extensions. It would appear on the left side of the workspace. 

When you complete the commit, you can monitor the test under _Actions_ bar at the top of github repository. If everything is working fine, then it will be finished without any errors! 


> **FINAL NOTE**: https://github.com/kubranarci/nf-core-mydemo/tree/dev includes a run-ready pipeline with a results directory. If you haven't managed to complete your workflow just yet, you can have a look. 

---


### Are you looking for your next challange? 

In this small exercise you will be creating a simple workflow to subset, sort and index given alignment files.

1. Create a new pipeline using _create_ command

```bash
nf-core create -n sortandindex -d "This pipeline subsets a region, sorts and indexes given alignment files" --plain
```

2. Edit schema_input.json to input **a BAM file** and **a string** to subset a region.

4. Set input_ch.

5. Create input channels for fasta/fai references. You can use igenomes in that purpose.

6. Install relevant samtools modules:

```bash
nf-core modules install samtools/sort
nf-core modules install samtools/view
nf-core modules install samtools/index

```
7. Include samtools modules to the workflow.

8. Connect modules so that an example run will perform the following tasks:

```
samtools view -o input.inputregion.bam -b input.bam inputregion

samtools sort -o input.inputregion.sorted.bam input.chinputregionr21.bam

samtools index -b  input.inputregion.bam input.inputregion.bam.bai

```

9. Edit modules.config to arrange arguments and output files properly

10. Edit test.config for quick test development. Tip: You can use the example bam files in created in this exercise.

11. Test your development

12. Once it is ready, commit your changes to your github repository.


### Sources

Who is GHGA?
- https://www.ghga.de/

Documentation and reference material for nextflow:
- Nextflow homepage: https://www.nextflow.io/
- Nextflow training material: https://training.nextflow.io
- Pipeline examples: https://www.nextflow.io/example1.html
- Main documentation: https://www.nextflow.io/docs/latest/index.html
- Common implementation patterns for developers: http://nextflow-io.github.io/patterns/index.html

Docker for beginners:
- https://docker-curriculum.com/
  
Github Actions & GitHub CI:
- https://docs.github.com/en/actions/automating-builds-and-tests/about-continuous-integration
- https://resources.github.com/ci-cd/


Do you have any questions? mail to kuebra.narci@dkfz-heidelberg.de 

