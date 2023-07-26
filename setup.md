---
title: "Data & Setup"
number-sections: false
---

<!-- 
Note for Training Developers:
We provide instructions for commonly-used software as commented sections below.
Uncomment the sections relevant for your materials, and add additional instructions where needed (e.g. specific packages used).
Note that we use tabsets to provide instructions for all three major operating systems.
-->

::: {.callout-tip level=2}
## Workshop Attendees

If you are attending one of our workshops, we will provide a training environment with all of the required software and data.  
If you want to setup your own computer to run the analysis demonstrated on this course, you can follow the instructions below.
:::

## Data

The data used in these materials is provided as a zip file. 
Download and unzip the folder to your Desktop to follow along with the materials.

<!-- Note for Training Developers: add the link to 'href' -->
<a href="https://www.dropbox.com/sh/0obk40tzxscqdez/AADJVWnZd3h8UoHlbYRq64wVa?dl=1">
  <button class="btn"><i class="fa fa-download"></i> Download</button>
</a>

## Setup

To run the analysis covered in this workshop, you will broadly need two things: 

- **R/RStudio** for all the downstream analysis (i.e. after peak calling using the `nf-core/chipseq` workflow). 
  These analyses can typically be run on your local computer and on any OS (macOS, Windows, Linux).
- A **Linux environment** to run the pre-processing steps and peak calling (i.e. running the `nf-core/chipseq` workflow). 
  We highly recommend using a dedicated server (typically a HPC) for this step. 
  Technically, you can also run this workflow on Windows via _WSL2_ (we provide instructions below), but we do not recommend it for production runs. 


### R and RStudio

::: {.panel-tabset group="os"}

#### Windows

Download and install all these using default options:

- [R](https://cran.r-project.org/bin/windows/base/release.html)
- [RTools](https://cran.r-project.org/bin/windows/Rtools/)
- [RStudio](https://www.rstudio.com/products/rstudio/download/#download)

#### macOS

Download and install all these using default options:

- [R](https://cran.r-project.org/bin/macosx/)
- [RStudio](https://www.rstudio.com/products/rstudio/download/#download)

#### Linux

- Go to the [R installation](https://cran.r-project.org/bin/linux/) folder and look at the instructions for your distribution.
- Download the [RStudio](https://www.rstudio.com/products/rstudio/download/#download) installer for your distribution and install it using your package manager.

:::


### R Packages

Open RStudio and run the following: 

```r
# install BiocManager if not installed already
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

# Install all packages used
BiocManager::install(c("GenomicRanges", "rtracklayer", "plyranges", "ChIPseeker", "profileplyr", "ggplot2", "DiffBind"))
```


### Conda/Mamba

For the command-line tools covered in the course you will need a Linux machine (or _WSL2_, if you are on Windows - see @sec-wsl).

If you are an experienced Linux user, you can install/compile each tool individually using your preferred method. 
Otherwise, we recommend doing it via the [Mamba package manager](https://mamba.readthedocs.io/en/latest/installation.html). 
If you already use Conda/Mamba you can skip this step. 

To make a fresh install of Mamba, you can run: 

```bash
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
```

And follow the instructions on the terminal, accepting the defaults. 
Make sure to restart your terminal after the installation completes.

These instructions also work if you're using a HPC server.


### Nextflow

We recommend having a dedicated environment for _Nextflow_, which you can use across multiple pipelines you use in the future. 
Assuming you've already installed Conda/Mamba, open your terminal and run: 

```bash
mamba create --name nextflow nextflow
```

Whenever you want to use nextflow, you need to activate your environment with `conda activate nextflow`.


### ChIP-seq tools

For other command-line tools that we covered in the workshop, you can install them in their own conda environment:

```bash
mamba create --name chipseq
mamba install --name chipseq idr deeptools meme homer
```

When you want to use any of them, make sure to activate your environment first with `conda activate chipseq`.


### Windows WSL2 {#sec-wsl}

:::{.callout-warning}
We highly recommend running the raw data processing pipeline on a dedicated Linux server (typically a HPC), not directly on Windows via WSL2.
Although you can technically run the entire pipeline on WSL2, it may be a very suboptimal way of doing so for real data.
:::

The **Windows Subsystem for Linux (WSL2)** runs a compiled version of Ubuntu natively on Windows. 
There are detailed instructions on how to install WSL on the [Microsoft documentation page](https://learn.microsoft.com/en-us/windows/wsl/install). 
Briefly:

- Click the Windows key and search for  _Windows PowerShell_, open it and run the command: `wsl --install`.
- Restart your computer. 
- Click the Windows key and search for _Ubuntu_, which should open a new terminal. 
- Follow the instructions to create a username and password (you can use the same username and password that you have on Windows, or a different one - it's your choice). 
- You should now have access to a Ubuntu Linux terminal. 
  This (mostly) behaves like a regular Ubuntu terminal, and you can install apps using the `sudo apt install` command as usual. 

#### Setup directories

After WSL is installed, it is useful to create shortcuts to your files on Windows. 
Your `C:\` drive is located in `/mnt/c/` (equally, other drives will be available based on their letter). 
For example, your desktop will be located in: `/mnt/c/Users/<WINDOWS USERNAME>/Desktop/`. 
It may be convenient to set shortcuts to commonly-used directories, which you can do using _symbolic links_, for example: 

- **Documents:** `ln -s /mnt/c/Users/<WINDOWS USERNAME>/Documents/ ~/Documents`
  - If you use OneDrive to save your documents, use: `ln -s /mnt/c/Users/<WINDOWS USERNAME>/OneDrive/Documents/ ~/Documents`
- **Desktop:** `ln -s /mnt/c/Users/<WINDOWS USERNAME>/Desktop/ ~/Desktop`
- **Downloads**: `ln -s /mnt/c/Users/<WINDOWS USERNAME>/Downloads/ ~/Downloads`


#### Docker for WSL

We've experienced issues in the past when running _Nextflow_ pipelines from WSL2 with `-profile singularity`. 
As an alternative, you can instead use **_Docker_**, which is another software containerisation solution. 
To set this up, you can follow the instructions given on the Microsoft Documentation: [Get started with Docker remote containers on WSL 2](https://learn.microsoft.com/en-us/windows/wsl/tutorials/wsl-containers).

Once you have _Docker_ set and installed, you can then use `-profile docker` when running your _Nextflow_ command.


### Singularity

Singularity is a software for running a virtual operating system locally (known as a container) and popularly used for complex bioinformatic pipelines. 
_Nextflow_ supports the use of _Singularity_ for managing its software and we **recommend its use it on HPC servers**. 
_Singularity_ is typically installed by your HPC admins, otherwise request that they do so. 

However, if you want to run the analysis locally on your computer (again, we do not recommend you to do so), then you can install Singularity following the instructions below.

::: {.panel-tabset group="os"}

#### Windows

You can use _Singularity_ from the _Windows Subsystem for Linux_ (see @sec-wsl).
Once you setup WSL, you can follow the instructions for Linux.

#### macOS

Singularity is [not available for macOS](https://docs.sylabs.io/guides/3.0/user-guide/installation.html#install-on-windows-or-mac).

#### Linux

These instructions are for _Ubuntu_ or _Debian_-based distributions[^1].

[^1]: See the [Singularity documentation page](https://docs.sylabs.io/guides/3.0/user-guide/installation.html#install-on-linux) for other distributions.

```bash
sudo apt update && sudo apt upgrade && sudo apt install runc
codename=$(lsb_release -c | sed 's/Codename:\t//')
wget -O singularity.deb https://github.com/sylabs/singularity/releases/download/v3.10.2/singularity-ce_3.11.4-${codename}_amd64.deb
sudo dpkg -i singularity.deb
rm singularity.deb
```

:::



<!-- 
### Visual Studio Code

::: {.panel-tabset group="os"}

#### Windows

- Go to the [Visual Studio Code download page](https://code.visualstudio.com/Download) and download the installer for your operating system. 
  Double-click the downloaded file to install the software, accepting all the default options. 
- After completing the installation, go to your Windows Menu, search for "Visual Studio Code" and launch the application. 
- Go to "_File > Preferences > Settings_", then select "_Text Editor > Files_" on the drop-down menu on the left. Scroll down to the section named "_EOL_" and choose "_\\n_" (this will ensure that the files you edit on Windows are compatible with the Linux operating system).

#### Mac OS

- Go to the [Visual Studio Code download page](https://code.visualstudio.com/Download) and download the installer for Mac.
- Go to the Downloads folder and double-click the file you just downloaded to extract the application. Drag-and-drop the "Visual Studio Code" file to your "Applications" folder. 
- You can now open the installed application to check that it was installed successfully (the first time you launch the application you will get a warning that this is an application downloaded from the internet - you can go ahead and click "Open").

#### Linux (Ubuntu)

- Go to the [Visual Studio Code download page](https://code.visualstudio.com/Download) and download the installer for your Linux distribution. Install the package using your system's installer.

:::
 -->
