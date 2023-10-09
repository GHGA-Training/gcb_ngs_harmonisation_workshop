## **Standardizing and harmonizing NGS analysis workflows**

Materials for NGS Harmonization Workshop at GCB 2023

It is planned to be 3 hours of an overview of standardization and harmonizing NGS analysis strategies in GHGA. We will explore how FAIR principles enable the standardization and harmonization of nf-core-based NGS analysis workflows within GHGA. We will  demonstrate the adaptability of nf-core workflows and discuss the importance of standardization of workflows. Finally, we will show how to make workflows scalable, robust, and automated using a small subset of a public dataset. 

### Preliminary Schedule

|Time|Topic|
|:---|:---|
|9:00am  - 9:15am|Introduction to the tutorial: What is GHGA? What are our workflow objectives? What is FAIR data|
|9:15am  - 9:45am|Reproducibility, adaptability, and portability of Workflows|
|9:45am  - 10:15am|Standardization of workflows using Workflow Managers|
|10:15am - 10:30am|Break|
|10:30am - 11:00am|Accurate analysis and benchmarking|
|11:00am - 12:00am|Hands-on experience (with minimal test cases provided)


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
- Raw files (one of them is enough):
    - [NA12878 chr21](https://drive.google.com/drive/folders/1OXGIx9RHioH1QB65SK75m_liP_fygxYH?usp=drive_link)
    - [SEQC2 small set](https://github.com/GHGA-Training/gcb_ngs_harmonisation_workshop/tree/main/reads) 

## Construction of a simple alignment and variant calling pipeline using nf-core tools

Please follow all of the steps 1-9  to be able to construct a simple nextflow pipeline. 


### 1. What is **GitHub**, and how we can use it?

The most common way of providing stable versions for software is through git systems. Github, gitlab or bitbucket are all valuable services in software development. They provide code hosting/storage platform for version control and collaboration. It lets you and others work on the same project from anywhere. Using git systems one can create different branches (or copies) of the stable working version dependent or independent of one another. 

GitHub is one of the most commonly used platform for collaborative software development and version control using git. We will be using the following repositories through this workshop. If you want, you can fork gcb_ngs_harmonisation_workshop or follow from the original documentation. _nf-core/testpipeline_ will be used as a template for our workflow construction. Again,  you can either fork and create a local version of it. Dont forget to create your working (dev) branch to work on. 

- fork https://github.com/GHGA-Training/gcb_ngs_harmonisation_workshop
- fork https://github.com/nf-core/testpipeline -r v0.1.5

Whenever you fork those pages, you can take any actions on that as you wish! 

```
git clone https://github.com/<yourgithubname>/testpipeline -r v0.1.5
git checkout -b dev
```

if you don't want to fork the repo, pull it in your local computer work on it and push it into your own github repo! 

### 2.  Creating your own VS-Code workstation:

There are plenty of good editors but here we will be using VS-Code. If you have any other preference just go with it! The idea is to be able to connect multi-services in an environment creating a functional development strategy for lining, debugging,  editing, and pushing your code. You can add, VS-code extensions like Nextflow, nf-core tools, github, python, R and many more. 

Simply, open _testpipeline_ using VS-code and start working on it!
  
### 3.  How to build and use **Docker** containers:

- Containers are configurable virtualization technology that allows the packaging and distribution of pipelines  in a self-contained and platform-independent manner. By this way, installed software using package managers can be packed into an image with its corresponding dependencies. An image can be considered as  a file (or set of files) that contains the application with the code and all its dependencies, libraries, etc. You can copy images around, upload them into registries, download them, and re-use them. Commonly used software for containerization in bioinformatics includes Docker, Singularity and podman. 

- Using software containers is crucial to make our pipelines portable. In this workshop, we will be using docker and singularity in order to dive in better, we will be constructing our own bcftools container as an example.  
  
- Download and install Docker [here]([url](https://docs.docker.com/get-docker/))

```
docker login
```

- The pull command lets you download a Docker image without running it. For example

```
docker pull ubuntu:focal
```

 - Launching a BASH shell in the container allows you to operate in an interactive mode in the containerized operating system. For example

```
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
```
docker commit <containerid> <docker_name>/imagename:version_tag
```

- Push the image to Dockerhub registry when the commit is done:

```
docker push <docker_name>/imagename:version_tag
```

- Singularity is another containerization platform that you can use instead of docker. Docker images can quickly be converted into singularity. It is also possible to use docker-converted singularity images directly by nextflow as follows:

singularity_image = docker://<docker_name>/imagename:version_tag


### 4.  Getting familiar with **Nextflow** and **nf-core**

The new generation way of pipelines are  Workflow management systems providing all components of standard analysis. They increase transparency, enable long-term sustainability, and aid in achieving findable, accessible, interoperable, and reusable computational analysis and they are FAIR-approved! Snakemake, Galaxy, cornwell, pegasus and Nextflow are all open source and free-to-use examples! Workflow managers use containers to run software providing portability and stability. They provide scalability: They automatically manage the scheduling enabling the effective utilization of available resources. Required reference files can be pulled through cloud registries ensuring saving space locally. They are able to cache intermediate results avoiding recalculations saving significant time and computing resources which is a key advantage of workflow managers. Project repositories are often in git systems allowing collaboration and feedback. When the execution is finished, They provide the option to generate execution reports with detailed information about the run. 

Nextflow is growing fast and has long-term support available from Seqera Labs. It has been developed since 2013 by the same team. Nf-core is a community effort to collect a curated set of analysis pipelines built using Nextflow. The effort itself is not just restricted to pipelines, they also provide sub-workflows and modules (single processes running tools). Their modular design is very standard and also flexible for multi-purposing.  They also provide helper tools to create template pipelines, install modules, and create test cases. 


```
nextflow info
```

**nf-core** is a community effort to optimize and collect a set of nextflow pipelines. Currently, they have **86** pipelines!

- let's explore **nf-core** tool:

```
nf-core --help
```

- Listing available pipelines:

```
nf-core list
```

- Listing available modules through nf-core:

```
nf-core modules list remote
```

### 5.  Using nf-core modules to build a simple variant calling pipeline

- Let's use **nf-core/testpipeline** as a template. We will use the local fork of the pipeline.
  
Note: Do not use _master_ brach switch into _dev_!

- We will perform _bwa-mem_ alignment, which requires indexed fasta genome, using _bwa-index_. Thus we need bwa-mem and bwa-index modules. Luckily, nf-core provides bwa modules and we can directly install them!

```
nf-core modules install bwa/mem
nf-core modules install bwa/index
```

Now both modules should be located in **modules/nf-core** directory. 

- In order to use those modules, we will add descriptions to the workflow. 

Add two lines of code to workflows/testpipeline.nf

```Nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BWA_MEM                } from '../modules/nf-core/bwa/mem/main'
include { BWA_INDEX              } from '../modules/nf-core/bwa/index/main'
```

- In order to perform alignment using _bwa-mem,_ we need the reference fasta file to be indexed. testpipeline includes igenome.config template ready to use. Therefore, we can directly use one of the provided fasta files readily. We will use _bwa_index_ tool to index fasta file. 

Note: [IGenomes](https://emea.support.illumina.com/sequencing/sequencing_software/igenome.html) is a source providing collections of references and annotations supported by AWS.

igenome.config includes parameters for the available sources for the pipeline. Yet, to be able to use them, we need to create a channel for them. Let's create a genome channel for the usage of BWA_INDEX module. 

```Nextflow
ch_genome_fasta = Channel.fromPath(params.fasta).map { it -> [[id:it[0].simpleName], it] }.collect()
```

- Now, we would only need to add our first module! Checking out BWA_INDEX, the only input file is a fasta file to be able to create bwa index directory.

The output index directory will be saved into ch_index channel for our further usage. 

```Nextflow
    //
    // MODULE: BWA_INDEX
    //
    BWA_INDEX(
        ch_genome_fasta
    )
    ch_index = BWA_INDEX.out.index
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
```

- Next task will be adding BWA_MEM alignment module:
- 
BWA_MEM module requires fastq reads which were already prepared through INPUT_CHECK. We prepared ch_index channel in the previous step.
"_true_" statement is a value in order to activate samtools sort for BAM file.

```Nextflow
    //
    // MODULE: BWA_MEM
    //
    BWA_MEM(
        INPUT_CHECK.out.reads,
        ch_index,
        "true"
    )
    ch_bam = BWA_MEM.out.bam
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)
```

- Final task will be calling variants using the bam files generated with BWA_MEM. We will use bcftools mpileup together with bcftools call and bcftools view. Let's create our first module using bcftools container that we created!

A draft BCFTOOLS_MPILEUP module will look like this: We will need to define input and output files together with bcftools commands that will process variant calling. 

``` Nextflow
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

```Nextflow
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

```Nextflow
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


```Nextflow
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
```

 We can add either our docker container (created in the 3rt step) or search for the proper bcftool version on[ biocontainers registry](https://quay.io/repository/biocontainers/bcftools?tab=tags). 

```Nextflow
    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'quay.io/biocontainers/bcftools:1.17--haef29d1_0' }"
```

 The final module should be like this:


```Nextflow
process BCFTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'quay.io/biocontainers/bcftools:1.17--haef29d1_0' }"

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

- Since we finished BCFTOOLS_MPILEUP module, we can include it in the workflow. First, add path of the module:

``` Nextflow
//
// MODULE: Installed directly locally
//
include { BCFTOOLS_MPILEUP            } from '../modules/local/bcftools_mpileup.nf'

```
  
- Let's connect BCFTOOLS_MPILEUP module to the output of BWA_MEM module

``` Nextflow
    //
    // MODULE: BCFTOOLS_MPILEUP
    //
    BCFTOOLS_MPILEUP(
        ch_bam,
        ch_genome_fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)
```

- Now, let's prepare our **input file** which includes our sample reads:

  - This pipeline comes with a ready sub-workflow in order to check the input files and create an input channel with them (subworkflows/input_check.nf). It automatically checks and validates the header, sample names, and sample directories. This module is implemented for both pair-end and single-end fastq file processing simultaneously. Since we will also use the same format, we won't change the module and make use of it directly. But, still, we need to prepare our own samplesheet accordingly to be able to input our files:

Place fastq files (reads) into testpipeline directory and create mysamplesheet.csv file: 


```console
sample,fastq_1,fastq_2
sample_paired_end,reads/NA12878_75M_Agilent_1.merged.fastq.gz,reads/NA12878_75M_Agilent_2.merged.fastq.gz
sample_single_end,reads/NA12878_75M_Agilent_1.merged.fastq.gz,
```


- Our simple pipeline, providing parallel alignments for both paired-end and single-end fastq files is ready! Now, we need to create a config file to describe the parameters needed for the run.

The config file needs to include minimal information about the run. 

Open a text file and create a config file as follows and name to mytest.config


``` Nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow config file for running minimal tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Defines input files and everything required to run a fast and simple pipeline test.

    Use as follows:
        nextflow run nf-core/testpipeline -profile mytest,<docker/singularity> --outdir <OUTDIR>

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
    // TODO nf-core: Specify the paths to your test data on nf-core/test-datasets
    // TODO nf-core: Give any required params for the test so that command line flags are not needed
    input  = 'assets/mysamplesheet.csv' 

    // Genome references
    genome = 'R64-1-1'
}
```

Then, we should add the path of our config ('conf/mytest.config' ) into nextflow.config file to be able to use directly. It should be inside _profiles_ section : 

``` Nextflow
profiles{
    mytest    { includeConfig 'conf/mytest.config'    }
}
```

One of the other configurations that we need to be careful of is setting a mirror to your container settings. In _nf-core/testpipeline_ they set _quay.io_ as a mirror for the registries by default. But if you want to use your own docker container generated in the 3rd step you should delete the default settings. 

``` Nextflow
apptainer.registry   = 'quay.io'
docker.registry      = 'quay.io'
podman.registry      = 'quay.io'
singularity.registry = 'quay.io'
```
to 

``` Nextflow
apptainer.registry   = ''
docker.registry      = ''
podman.registry      = ''
singularity.registry = ''
```

If you decide to use ready containers from biocontainers default settings can remain. 

### 6.  Our first full-functioning pipeline is ready! and we can directly run it!

Here, we will try to run our workflow using the config file we created (mytest), and input file we prepared (mysamplesheet.csv) and docker environment. But you can also run it with singularity. You can define another output directory else than _results_. 
```
nextflow run main.nf -profile mytest,docker --outdir results --input mysamplesheet.csv
```

We can actually test and debug our pipeline using this command. What is really cool and helpful to use the **-resume** tag to resume previously finished jobs! Moreover, don't forget to check out .nextflow.log files in case of an error. All of the runs will be saved into _work_ directory. 

In the first run of this workflow, container images and genome would need to be pulled (all automatic) which is why the initial run will take a bit longer than the following runs!  

NOTE: If you receive a permission denied error for the bin directory just provide the necessary permissions as follows:
```
chmod +x bin/*
```

### 7. Analyzing the results:

- _Provenance_ describes the history of a computational experiment and it is a very important aspect of workflows. Minor changes, usually in the workflow environment, software versions, parameter settings, and reference annotations, may or may not change the results of each workflow run but the major challenge is to know and record the exact environment you run your analysis. Data provenance describes this trail of methods, versions, and arguments that were used to generate a set of files. In short, it is the history of a computational environment you used for the analysis.
- nf-core pipelines are armed through many provenance tools with the construction of _work_ directory with _temp_ results of each run, saving workflow execution history through .nextflow.log files as well as _CUSTOM_DUMPSOFTWAREVERSIONS_ module and MultiQC tools. 
- Collection of versions is a vital process in order to keep track of software history. In nf-core pipelines, each tool version is collected in a channel and then processed using _CUSTOM_DUMPSOFTWAREVERSIONS_ module and represented through MultiQC tool.
- Used data, produced inputs and outputs, artifacts, information about the execution like consumed memory and CPU, processed time saved and represented in final reports of nf-core tools. 
- MultiQC tool also aggregates logs and reports from the analysis. In our analysis, FASTQC analysis was already included. In this example file, you can both see the FASTQC report and the software versions together with the workflow summary.


### 8. Committing changes to github repo, and Continuous Integration test environment  

Another advantage of using github is easy debugging. Testing the code and building can be automized. Continuous Integration (CI) tests refer to the build and unit testing stages of the software release process and are curious for automatic software development. Every revision that is committed triggers an automated build and test. When it is reviewed and approved, it can be merged safely. The aim is to protect the main branch from errors and unintended changes. CI tests are mandatory for nf-core pipelines and luckily in our example nf-core/testpipeline also comes with an environment that we can actually work on! In order to activate automatic tests:

- Activate actions tab on github (on your github repo testpipeline)
- under .github/workflows/ci.yml

Change this line

```yml
    # Only run on push if this is the nf-core dev branch (merged PRs)
    if: "${{ github.event_name != 'push' || (github.event_name == 'push' && github.repository == 'nf-core/testpipeline') }}"
```
to
```yml
    # Only run on push
    if: "${{ github.event_name != 'push' || (github.event_name == 'push' && github.repository == '<your github name>/testpipeline') }}"
```

Now whenever we commit changes to _dev_ branch, a small test pipeline will be triggered in order to check if the workflow/code works smoothly.  

    
The last thing we need to do is to commit the changes we have in our repository and push into _dev_ branch that we were working in. 

```
git add .
git commit -a -m "commit" (commit message)
```

It is also possible to commit changes using VS-Code if you installed github extensions. It would appear on the left side of the workspace. 

When you complete the commit, you can monitor the test under _Actions_ bar at the top of github repository. If everything is working fine, then it will be finished without any errors! 


**FINAL NOTE**: https://github.com/kubranarci/testpipeline/tree/dev includes a run-ready pipeline with a results directory. If you haven't managed to complete your workflow just yet, you can have a look. 


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

