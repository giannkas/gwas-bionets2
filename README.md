# gwas-bionets2
Network-based approaches for gene discovery by leveraging GWAS data

In comparison with previous version [gwas-bionets](https://github.com/giannkas/gwas-bionets.git), 
this version is developed using Nextflow domain-specific language 2, and includes more
network-based methods: Hierarchical HotNet, SigMod, Heinz, dmGWAS, LEAN and HotNet2 (you may exclude it if Hierarchical HotNet is used).
Although gwas-bionets resembles a previous version of this repository, I decided not to continue in the same repo simply because the technologies are different. I want to ensure that everything is well-documented in case someone needs to reproduce it, especially in HPC environments where permissions are restricted and programs are often outdated by default. Also, I have not decided to include everything in a container because not all of the methods can be freely distributed to third parties.

## Software requirements

I explain how to install all the network-based methods, and how to set the environment for running them.

Ideally, create a folder in your home directory to store all software. For example:

```bash
mkdir ~/bin
```

**Install Java (required for installing nextflow)**

Create a "java" folder in the software directory and navigate to it.

```bash
mkdir ~/bin/java
cd ~/bin/java
```

Download the `x64` version of Java JDK as a `tar.gz` file from the Oracle website into your machine and decompress it. In this case is version 24, but you can check [latest releases](https://www.oracle.com/java/technologies/downloads/).


```bash
wget https://download.oracle.com/java/24/latest/jdk-24_linux-x64_bin.tar.gz
tar xzfv jdk-24_linux-x64_bin.tar.gz 
rm jdk-24_linux-x64_bin.tar.gz
```

Export the path to the `bin` directory of this folder into the system variable `$PATH` to make Java executable. Also, export the `$JAVA_HOME` variable indicating the root directory. Ideally, add these to `~/.bashrc` to avoid repeating the process on each server connection or reboot, eg.

```bash
export PATH=/home/username/bin/java/jdk-24.0.1/bin:$PATH
export JAVA_HOME=/home/username/bin/java/jdk-24.0.1
```

You may need to source `.bashrc` file before checking installation, so type:

```bash
source ~/.bashrc
```

Test the installation:

```bash
java -version
```

You should see something like:

```bash
java version "24.0.1" 2025-04-15
Java(TM) SE Runtime Environment (build 24.0.1+9-30)
Java HotSpot(TM) 64-Bit Server VM (build 24.0.1+9-30, mixed mode, sharing)
``` 

**Install Nextflow**

People at [sequera](https://seqera.io/) provide better explanations and documentation about Nextflow than I could offer. Please follow the [installation steps there](https://www.nextflow.io/docs/latest/install.html#install-nextflow).

**Install MAGMA**

You can follow the installation instructions for MAGMA at its website (version 1.10): [Multi-marker Analysis of GenoMic Annotation](https://cncr.nl/research/magma/). According to the documentation MAGMA: is a self-contained executable and does not need to be installed. 

Go to our software directory and create a 'magma' folder where the binaries will be located.

```bash
cd ~/bin
mkdir magma
```

Download the `zip` file.

```bash
curl -v "https://vu.data.surfsara.nl/index.php/s/zkKbNeNOZAhFXZB/download" -H "Accept-Encoding: zip" > magma_v1.10.zip
```

Descompress `magma_v1.10.zip` file. We will decompress to the created folder `magma` (unsing `unzip -d` option)

```bash
unzip magma_v1.10.zip -d magma
```

Remove the `zip` file.

```bash
rm magma_v1.10.zip
```

Test MAGMA

```bash
./magma/magma
```

You should see something like:

```bash
No arguments specified. Please consult manual for usage instructions.

Exiting MAGMA. Goodbye.
```

You can decide to include MAGMA's location into the PATH variable so it is called system-wide under your session. Otherwise, you must indicate where to find MAGMA when calling `magma_calc.nf` script; it has a parameter named `magma` for this purpose.

**Install PLINK**

Similarly, you can install PLINK (version 1.9) from its website: [population linkage](https://www.cog-genomics.org/plink/1.9/). PLINK is also self-contained executable so either you add to your path or reference the executable when using it.

Again, please locate our software directory and create a 'plink' folder where the binaries will be saved.

```bash
cd ~/bin
mkdir plink
```

Download the `zip` file

```bash
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip
```

Decompress `plink_linux_x86_64_20241022.zip` file. We will decompress to the created folder `plink` (unsing `unzip -d` option)

```bash
unzip plink_linux_x86_64_20241022.zip -d plink
```

Remove the `zip` file.

```bash
rm plink_linux_x86_64_20241022.zip
```

Test PLINK

```bash
./plink/plink --version
```

You should see something like:

```bash
PLINK v1.9.0-b.7.7 64-bit (22 Oct 2024)
```

We have to include PLINK's location into our PATH variable because there is no parameter in the pipeline to reference its location. Conversely, we did include a PLINK parameter to indicate PLINK version, either 1 or 2. Then, to add plink to the environment variable, you proceed as follow:

```bash
export PATH=$PATH:/home/username/bin/plink
```

You can add it to your `.bashrc` file to make it permanent. And then, `source ~/.bashrc` to apply changes in your current session.

**Install R and some packages (required for the methods)**

If you dont't have R installed in your machine (add your superuse credentials if needed to install software), then proceed as follows (or you can check instructions from [The Comprehensive R archive Network](https://cloud.r-project.org/)):

```bash
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
apt-get -y install --no-install-recommends r-base r-base-dev
```

Setup the general CRAN repo.

```bash
echo 'local({
    r <- getOption("repos")
    r["CRAN"] = "https://cloud.r-project.org/"
    options(repos = r)
  })' >> /etc/R/Rprofile.site
```

Install Bioconductor (`BiocManager`), `twilight` and `BioNet` (the latter contains the necessary files to use Heinz method).

```bash
R -e "install.packages('BiocManager')"
R -e "BiocManager::install('BioNet')"
R -e "BiocManager::install('twilight')"
```

Install R packages, `tidyverse`, `cowplot`, `igraph` and `gprofiler2`:

```bash
R -e "install.packages(c('tidyverse', 'cowplot', 'igraph', 'gprofiler2', 'foreach'))" 
```

**Install LEAN**

As of the time of writing this README, LEANR has been [removed](https://cran.r-project.org/web/packages/LEANR/index.html) from the CRAN repository. Therefore, a manual installation is recommended.

```bash
wget https://cran.r-project.org/src/contrib/Archive/LEANR/LEANR_1.4.9.tar.gz
```

And then install it using `R console`.

```R
install.packages("path/to/LEANR_1.4.9.tar.gz", repos = NULL, type = "source")
```

**Install dmGWAS**

The pipeline uses dmGWAS v3.0 or EW_dmGWAS released on October 4, 2014, one can refer to this [page](https://bioinfo.uth.edu/dmGWAS/) for more information. Download it as follows:

```bash
wget https://bioinfo.uth.edu/dmGWAS/dmGWAS_3.0.tar.gz
```

And then install it using `R console`

```bash
install.packages("path/to/dmGWAS_3.0.tar.gz", repos = NULL, type = "source")
```

**Install SigMod**

You can install SigMod (version 2) from this website: [Strongly Interconnected Gene MODule](https://github.com/YuanlongLiu/SigMod/tree/20c561876d87a0faca632a6b93882fcffd719b17). It contains R scripts, then it suffices to assign the parameter `sigmod_path` when calling the `bionets.nf` script, eg. `--sigmod_path="~/bin/SigMod_v2"`.

You change directory to our software directory `bin` folder:

```bash
cd ~/bin
```

Download the `zip` file and decompress it (no need to create a folder since there is a folder inside including the code and manual):

```bash
wget https://github.com/YuanlongLiu/SigMod/raw/20c561876d87a0faca632a6b93882fcffd719b17/SigMod_v2.zip
unzip SigMod_v2.zip
```

Change folder name

```bash
mv SigMod_v2 sigmod
```

**Install Heinz**

You have already installed it when installing the `BioNet` package from Bioconductor :-)

After that, we are all set!

## Main Scripts

1. This script works with the raw data for splitting it if parametrized with the _k_ parameter. The script has the needed parameters to be filled by the user, clearly, you can run each of the steps within the script separately.

`bionets_construction_from_data.sh`

2. This script works with the scores previously computed using a software like [MAGMA](https://cncr.nl/research/magma/) for the gene P-values. Again, a _k_ parameter greater than 1 generates k-fold solutions. As above, we conceived the script to be modified to provide the parameters.

`bionets_construction_from_scores.sh`
