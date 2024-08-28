# Exploring Analysis  Platforms and High-Performance Computing Resources for Bioinformatics 

<br>
<br>
<p align="center">
    <img src="https://hpcdocs.hpc.arizona.edu/assets/images/home/uofa_cactus_computing3.png" width="300">
</p>
<br>

---
>[!important]
> :clock1: **Schedule**
> - 2:00pm-2:10pm: Welcome and overview of the DataLab Biosciences workshop
> - 2:10pm-2:45pm: Introduction to the HPC
> - 2:45pm-3:00pm: Introduction to CyVerse

>[!important]
> :heavy_exclamation_mark: **Requirements**
> - Basic command line knowledge
>- Access to a [Terminal](https://en.wikipedia.org/wiki/Unix_shell)
>    - Unix and Mac users already have access to the Terminal
>    - Windows users can use either [PowerShell](https://en.wikipedia.org/wiki/PowerShell) or the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install)
> - A registered CyVerse account (Register for a CyVerse account)

>[!important]
> :white_check_mark: **Expected Outcomes**
> - Understanding the structure of an HPC
> - Learning basic HPC commands
> - Exposition to on-campus Cloud computing alternatives (CyVerse)

<br>

---
---

## Bioinformatics and the Need for Compute

The field of bioinformatics has become increasingly reliant on computational power. The volume and complexity of biological data (e.g., genome sequences, protein structures and molecular simulations) requires advanced computational resources for data extrapolation and analysis.

To effectively tackle the problem, scientists can access resources such as [High Performance Computing](https://en.wikipedia.org/wiki/High-performance_computing) (HPC) systems or [cloud based computing](https://en.wikipedia.org/wiki/Cloud_computing) to tackle the issue at hand.

The University of Arizona allows for access to both of these platforms through the U of A HPC system (https://hpc.arizona.edu/) or CyVerse (https://www.cyverse.org/). Both of these offer a number of software and resources that can help researchers in their scientific journey.

In this workshop, we are going to demonstrate the basics of connecting and using the resources here on Campus with the aim to tackle your computational needs. **Strap in!**

---
---

## A 10,000ft View of the HPC

<p align="center">
    <img src="https://hpcdocs.hpc.arizona.edu/quick_start/what_is_hpc/images/simple_hpc_diagram.png" width="450">
</p>

Here, we discuss how University of Arizona users can access the HPC and how to use the resources it offers. This section covers:
- [Logging onto the HPC](#logging-onto-the-hpc)
- [Choosing a system](#choosing-a-system) 
- [Checking Available Resources](#checking-available-resources)
- [Job Submissions](#job-submissions)

### Logging onto the HPC

If you have a UA account, to connect to the HPC you need to use `ssh` ([Secure Shell](https://en.wikipedia.org/wiki/Secure_Shell)). Open a terminal, and type:

```
ssh <UA username>@hpc.arizona.edu
```

Type your UA password and if successful you'll be greeted with a two-factor login. Select which choice, and complete the authentification. Once you are past the authentification steps, you will enter the [Bastion server](https://en.wikipedia.org/wiki/Bastion_host). This step has 2 purposes: 

1. Protect from attacks.
2. Select what HPC system you want to use.

> [!NOTE]
> The Bastion server is NOT YET the HPC! Here you cannot submit jobs or run analyes. Type `shell` in order to select what system you want to use.

The whole process (from logging to selecting the system) looks like the following:

```
ssh cosi@hpc.arizona.edu
(cosi@hpc.arizona.edu) Password: 
(cosi@hpc.arizona.edu) Duo two-factor login for cosi

Enter a passcode or select one of the following options:

 1. Duo Push to XXX-XXX-8418
 2. SMS passcodes to XXX-XXX-8418

Passcode or option (1-2): 1
Success. Logging you in...
Last login: Tue Mar 26 14:52:39 2024 from dhcp-10-132-212-1.uawifi.arizona.edu
This is a bastion host used to access the rest of the RT/HPC environment.

Type "shell" to access the job submission hosts for all environments
-----------------------------------------

[cosi@gatekeeper ~]$ shell
Last login: Wed Mar 20 10:30:25 2024 from gatekeeper.hpc.arizona.edu
***
The default cluster for job submission is Puma
***
Shortcut commands change the target cluster
-----------------------------------------
Puma:
$ puma
(puma) $
Ocelote:
$ ocelote
(ocelote) $
ElGato:
$ elgato
(elgato) $
-----------------------------------------

[cosi@wentletrap ~]$ ocelote
(ocelote) [cosi@wentletrap ~]$
```

At this point you are in the Login Node, where you can submit jobs or ask for an interactive node.

<p align="center">
    <img src="https://uarizona.atlassian.net/wiki/download/thumbnails/75989999/HPCDiagram_FileTransfers.png?version=1&modificationDate=1696282205000&cacheVersion=1&api=v2&effects=drop-shadow&width=1124&height=686" width="600">
</p>

### Choosing a System

In the example above, we chose the Ocelote system. Notice how there are 2 other choices: Puma and El Gato. Without going in  too much detail, here are some of the statistics regarding the 3 HPC systems at the UA HPC.

|System|Year of Aquisition|Processors|RAM|GPU|
|-|-|-|-|-|
|Puma|2020|2x AMD Zen2 48 CPU (94 cores total)|512GB|6x Nvidia V100S|
|Ocelote|2016|2x Xeon E5-2695v3 14-core (28 cores total)|192GB|46x Nvidia P100|
|El Gato|2013|2x Xeon E5-2650v2 8-core (16 core total)|64GB|removed as obsolete|

Find the full systems specs at the [official UA HPC documentatio resources page](https://uarizona.atlassian.net/wiki/spaces/UAHPC/pages/75990208/Compute+Resources).

El Gato is the oldest system, and potentially not useful for heavy research. Puma is the newest and most requested, whislt Ocelote is the "middle child": not as popular but still able to pack a punch. 

Depending on what your work is, your best bet would be Puma for heavy computation (if you are ok with waiting long queues); However, if your jobs aren't as taxing, then Ocelote could easily be a safe choice.

For this workshop, we are going to be using Ocelote.

### Checking Available Resources

#### Allocations

Once past the Bastion server and logged into the Ocelote login node, you are able to submit jobs or request an interactive node. Before you do so, it is wise to check your available resources. These resources are the ones made available to you by your PI or working group.

In order for you to check your resources, type `va`.

```
(ocelote) [cosi@wentletrap ~]$ va
Windfall: Unlimited

PI: parent_1743 Total time: 100000:00:00
	Total used*: 1:00:00
	Total encumbered: 0:00:00
	Total remaining: 99999:00:00
	Group: datalab Time used: 1:00:00 Time encumbered: 0:00:00


*Usage includes all subgroups, some of which may not be displayed here
```

`va` allows you to view all of the resources available from all the groups you are part of.

#### Storage

There are a number of ways one can approach storage on the HPC:

- Your own folder (in `/home/`): 50GB limit
- Your group (in `/groups/`): 500GB limit
- Your PI research (in `/xdisk/`): 20TB

Four the purpose of this workshop, we can access the datalab group storage in `/groups/cosi` (we will try and rename it in the future).

#### Queues and Submissions

One can check queue times and sumbissions by executing the SLURM command `squeue`. This will display the job id, submission type (standard, windfall, high priority), name of submission, user, status (queued, running), time elapsed, number of used nodes, and nodelist.

```
(ocelote) [cosi@wentletrap ~]$ squeue
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
2875825_[432-1000]  standard    qchem snanayak PD       0:00      1 (AssocGrpCPUMinutesLimit)
           2890485  standard R92L_Met krishnag PD       0:00      1 (Dependency)
           2890484  standard R92L_Met krishnag PD       0:00      1 (Dependency)
           2890483  standard R92L_Met krishnag PD       0:00      1 (Dependency)
           2881573  standard R92W_Met krishnag PD       0:00      1 (Dependency)
           2881572  standard R92W_Met krishnag PD       0:00      1 (Dependency)
           2802511  standard    eigen krishnag PD       0:00      1 (Dependency)
           2949224  standard brutefor theronir PD       0:00      1 (None)
           2898419  standard start_so mmartino PD       0:00      1 (Dependency)
           2898418  standard make_sor mmartino PD       0:00      1 (Dependency)
           2898416  standard start_so mmartino PD       0:00      1 (Dependency)
           2898415  standard make_sor mmartino PD       0:00      1 (Dependency)
           2898410  standard start_so mmartino PD       0:00      1 (Dependency)
           2898409  standard make_sor mmartino PD       0:00      1 (Dependency)
              .         .        .        .     .       .         .   .
              .         .        .        .     .       .         .   .
              .         .        .        .     .       .         .   .
           2884142  windfall li_d_t_6    bubin  R 2-06:58:50      1 i5n14
           2884107  windfall li_d_t_4    bubin  R 2-07:00:26      1 i6n8
           2884098  windfall li_d_t_7    bubin  R 2-07:00:50      1 i6n5
           2880486  windfall   be_10b   teodar  R 4-22:35:44      1 i16n5
           2880487  windfall    be_9b   teodar  R 4-22:35:44      1 i16n6
           2880488  windfall    be_7b   teodar  R 4-22:35:44      1 i16n12
           2880489  windfall    be_2b   teodar  R 4-22:35:44      1 i16n16
```

Likely, this will output 100s of lines, therefore if you want to check on your own job, you could use the CLI and `grep` to select the running submission (e.g., `squeue | grep <username>` or `squeue --user $NETID`).

---

## HPC, SLURM and Jobs Submissions

<p align="center">
    <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/3a/Slurm_logo.svg/1200px-Slurm_logo.svg.png" width="300">
</p>

All of the UA HPC systems run on a workload manager and job scheduler named [SLURM](https://slurm.schedmd.com/documentation.html) (Simple Linux Utility for Resource Management). It's designed to manage and schedule computing resources such as CPUs, GPUs, memory, and storage across a cluster of interconnected nodes.

You can learn more on SLURM and HPC system commands [here](https://uarizona.atlassian.net/wiki/spaces/UAHPC/pages/75989875/Running+Jobs+with+Slurm). 

### Job Submissions
There are 2 ways one can submit jobs onto the HPC system. The first is to run a **batch job**, which is the more popular submission type, whilst the other is by requesting an **interactive node**.

#### Batch jobs

As we are not going to be using batch submissions, we are not going to be going into too much detail. However, here is what you need to know. For more details on running batch jobs, visit the [official documentation page on batch jobs](https://uarizona.atlassian.net/wiki/spaces/UAHPC/pages/75989875/Running+Jobs+with+Slurm).

##### Writing a Batch Script

Batch scripts require a number of job **directives**. These are similar to the Dockerfile instructions, but instead of telling Docker how to build the image, these instead tell the SLURM system what to do with the job. The essential directives are the following:

|Directive|Purpose|
|-|-|
| `#SBATCH --account=group_name`	| Specify the account where hours are charged. |
| `#SBATCH --partition=partition_name` |	Set the job partition. This determines your job's priority and the hours charged. |
| `#SBATCH --time=DD-HH:MM:SS` |	Set the job's runtime limit in days, hours, minutes, and seconds. A single job cannot exceed 10 days or 240 hours. |
| `#SBATCH --nodes=N`	| Allocate N nodes to your job. |
| `#SBATCH --cpus-per-task=M` <br> and <br> `#SBATCH --ntasks=N` | ntasks specifies the number of tasks (or processes) the job will run. By default, you will be allocated one CPU/task. This can be increased by including the additional directive --cpus-per-task. |
| `#SBATCH --mem=Ngb`	| Select N gb of memory per node. If "gb" is not included, this value defaults to MB. |

After setting your directives, you can instruct the HPC to do what you require similar to a bash script.

Here's an example of a batch job:

```
#!/bin/bash
#SBATCH --job-name=blast_job          # Job name
#SBATCH --partition=standard          # Sets the job priority to standard
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks=1                    # Number of tasks (processes) per node
#SBATCH --cpus-per-task=4             # Number of CPU cores per task
#SBATCH --mem=8G                      # Memory per node (in this case, 8GB)
#SBATCH --time=02:00:00               # Time limit (HH:MM:SS)

# Load necessary modules
module load blast/2.12.0              # Load BLAST module (adjust version as needed)

# Change to the directory where the job will run
cd $SLURM_SUBMIT_DIR

# Define input and output files
query_file="query.fasta"              # Input query file (FASTA format)
database_file="database.fasta"        # BLAST database file (FASTA format)
output_file="blast_results.out"       # Output file for BLAST results

# Run BLAST command
blastp -query $query_file -db $database_file -out $output_file -evalue 0.001 -num_threads $SLURM_CPUS_PER_TASK
```

##### Submitting a Batch Script

- To submit jobs you need to use `sbatch`, such as `sbatch script.slurm`
- To cancel your job you do `scancel`, such as `scancel $JOBID` or `scancel -u $NETID`

This will submit your job to the queue. Execution will depend on your submission type (partition).

#### Launching an Interactive Node

An **interactive node**, unlike batch jobs which are run asynchronously, allows immediate access to compute. Similar to batch jobs, interactive nodes are submitted to the queue, but once available, you will receive a prompt for a node with the selected resources. Read more on how to launch interactive jobs in [the official documentation](https://uarizona.atlassian.net/wiki/spaces/UAHPC/pages/75989825/Interactive+Jobs).

> [!TIP] 
> **The "Quick and Dirty"**
>
> Don't need a lot of resources and just want access to the compute?
> Just type `interactive`. 
> Disclaimer: you may require to wait longer as your job is going to fall in the `windfall` queue.

Following are a list of useful flags (options) for setting up the interactive node.

|Flag|Default value| Description| Example|
|-|-|-|-|
| `-n` |	1 |	Number of CPUs requested per node	| interactive -n 8 |
| `-m` | 4GB |	Memory per CPU | interactive -m 5GB |
| `-a` |	none |	Account (group) to charge |	interactive -a datalab |
| `--partition=` |	windfall	| Partition to determine CPU time charges and is set to windfall when no account is specified, and is **set to standard when an account is provided.** | interactive --partition=windfall |
| `-t` |	01:00:00 |	Time allocated to session. | interactive -t 08:00:00 |
| `-N`	| 1 |	Number of nodes. | There is no reason to exceed 1 node unless the number of CPUs requested is greater than the number of CPUs per node on a given cluster.	| (in el gato, where number of CPU per node is 16 max ) interactive -n 32 -N 2 |

An example for an interactive node is:

```
interactive -n 8 -m 16GB -a datalab -t 02:00:00 
```

The above example will request an interactive node with 8 cores, 16GB RAM, "charging" the datalab, running for 2 hours. Try it!  

> [!TIP] 
> **Modules**
> There are 100s of tools installed on the HPC, few of which are available on the login screen. These tools are available only during a batch job submission or within interactive jobs.
>
>To see what tools are already running, or which are available, you will need to use the `module` command.

> [!TIP] 
> **Helpful `module` command**
> 
> |Command|Purpose|
> |-|-|
> |`module list`| Lists loaded modules|
> |`module avail`| Lists available modules|
> |`module spider`| Lists ALL modules|
> |`module load`| Loads a module|
> |`module help`| Help command!|

<br>

---
---

## CyVerse: Your Virtual Lab Bench

<br>
<p align="center">
    <img src="https://cyverse.org/sites/default/files/cyverse_logo_1_0.png" width="450">
</p>
<br>

[CyVerse](https://www.cyverse.org/) is an NSF funded cloud based system that aims to support scientists and researchers by allowing users to harness the power of cloud systems whilst making the cloud feel like your laptop.

Similar to the HPC, users can access a number of software useful for your own research and data extrapolation. This section covers the following content:
- [Accessing the CyVerse Discovery Environment](#accessing-the-cyverse-discovery-environment)
- [Data Management](#data-management-with-cyverse)
- [Running Applications and Resource Managenent](#running-applications-and-resource-managenent)
- [Sharing Data](#sharing-data-and-apps)

> [!WARNING] 
> As of 2022, CyVerse has implemented a [Subscription sytem](https://cyverse.org/subscribe); However, as part of the University of Arizona, students, faculty and staff can freely access CyVerse at the equivalent of the Pro-tier subscription

> [!WARNING]
> For your best experience, we suggest creating a CyVerse account at https://user.cyverse.org/ :warning:

### Accessing the CyVerse Discovery Environment

> [!TIP] 
> login to the Discovery Environemt at https://de.cyverse.org/ and click the top right Log In button!

One of CyVerse's main strengths is its [Discovery Environemt](https://de.cyverse.org/): a platform where users can launch Apps, view past and current Analyses, carry out Data Management and view avalable resources.

<p align="center">
    <img src="https://learning.cyverse.org/assets/de/de_login_screen.png" width="750">
</p>

The left hand side menu expands in the following:

- <img src="https://learning.cyverse.org/assets/de/menu_items/homeIcon.svg" width="15"> **Home/Dashboard**: Your main control panel that may display summary widgets, quick links to recent activities, or educational content such as tutorials and webinars.
- <img src="https://learning.cyverse.org/assets/de/menu_items/dataIcon.svg" width="15"> **Data**: This interface connects you to the Data Store. Here, you can manage your files, including uploading, downloading, organizing, and sharing data. You'll have access to your personal storage space and shared directories.
- <img src="https://learning.cyverse.org/assets/de/menu_items/appsIcon.svg" width="15"> **Apps**: Discover various applications, including VICE (Visual Interactive Computing Environment) apps for interactive computing sessions. You can browse, search, and launch these applications based on your research needs.
- <img src="https://learning.cyverse.org/assets/de/menu_items/analysisIcon.svg" width="15"> **Analyses**: View and manage your computational tasks. This section logs your history of analysis jobs, allowing you to monitor current processes, review completed ones, and access resulting data.
- <img src="https://learning.cyverse.org/assets/de/menu_items/webshellIcon.svg" width="15"> **Cloud Shell**: Access a Linux shell environment directly within the DE. This feature enables advanced users to perform command-line operations without leaving the platform.
- <img src="https://learning.cyverse.org/assets/de/menu_items/teamsIcon.svg" width="15"> **Teams**: Create and manage collaboration groups. Teams allow you to group together with other users for easier sharing of data, analyses, and other collaborative efforts.
- <img src="https://learning.cyverse.org/assets/de/menu_items/bank.svg" width="15"> **Collections**: Explore public collections of data and apps curated by other users or the CyVerse team. This resource can be invaluable for finding information relevant to your studies.
- <img src="https://learning.cyverse.org/assets/de/menu_items/helpIcon.svg" width="15"> **Help**: Access various support materials, including FAQs, guides, and contact information for direct assistance from the CyVerse support team.

The <img src="https://learning.cyverse.org/assets/de/menu_items/homeIcon.svg" width="15"> **Home/Dashboard** will also feature a quick breakdown of:
- Completed and recent Analyes
- Running Apps notifications
- Available resources (storage, compute time)

From here, users are able to manage data, launch Apps, and manage a variety of resources. 

### Data Management with CyVerse

Opening the <img src="https://learning.cyverse.org/assets/de/menu_items/dataIcon.svg" width="15"> **Data** tab, will show you your own private data. **This data is only visible to you and no one else**. [Sharing data](#sharing-data-and-apps) and Apps is covered below. 

Listed below are the methods to upload (and download) data with CyVerse:

| Method | Access Point	| Upload/Download | Installation/Setup Required | Account Required |	Max File Size |
|---|---|---|---|---|---|
| Discovery Environment	| Web |	Both | No |	Yes | 2GB/file upload, no limit import |
| Data Commons | Web | Download | No | No | 2GB/file |
| Cyberduck | Desktop App | Both | Yes | Yes, or public data only | >10GB |
| iCommands | Command line | Both | Yes | Yes, or public data only | >10GB |
| SFTP | Desktop App & Command line | Both | No (Yes for desktop apps) | Yes, or public data only | >10GB |
| GoCommands | Command line	| Both | Yes | Yes, or public data only	| >10GB |

Newcomers are encouraged to use any of the following methods:
- The Discovery Environement upload and download functions (useful for small files)
    - Note: the upload button is on the top right of the Data portal, whilst the download button will show only when a file/folder is selected
- [CyberDuck](https://cyberduck.io/): a great and easy way to access data in CyVerse as it allows for easy drag and drop from your own machine
- [GoCommands](https://learning.cyverse.org/ds/gocommands/): if you are versed with the Command Line, GoCommands are an in-house system for moving data as it is quicker than other methods and presents a number of additional options implementable in pipelines and automation. 

### Running Applications and Resource Managenent

Apps on CyVerse are accessible in the Apps tab on the left hand side menu. CyVerse supports a number of popular Apps such as Jupyter, RStudio, VSCode, Remote Desktops and a CLI. 

#### Launching and terminating an App

<p align="center">
    <img src="https://learning.cyverse.org/assets/de/de_apps_screen.png" width="750">
</p>

In order to launch an App, go to the App page in the Discovery Environment.

1. Select the App you want to launch by clicking on the name; This is going to open the App launch page/process.
2. Select the version, (if needed) the output folder, and give the Analysis a name and a comment. Click Next.
3. If the analysis is an executable, fill and give any required optional fields.
4. Click Submit. The page is going to turn from the App launch to Analysis.

From here on, you can access your analysis from the Analyes menu tab. Here you can go to the runnig App, share, or Terminate your analysis (by hitting the :x: button).

> [!TIP] 
> Terminating Analyses will save you compute hours!!

<p align="center">
    <img src="https://learning.cyverse.org/assets/de/menu_items/terminate_analyses.png" width="750">
</p>

#### Bring your own App

Users are welcome to bring their own Apps: Apps on CyVerse are [Containerized Images](https://www.docker.com/resources/what-container/) which have been implemented into CyVerse through an ingestion process:

1. User creates or obtains a container for a specific software. Containers are available at popular repositories like [DockerHub](https://hub.docker.com/) or [Quay.io](https://quay.io/).
2. Said Container is then made into a Tool on CyVerse, a process that allows to specify which [ports](https://en.wikipedia.org/wiki/Port_(computer_networking)) the software needs to be open, set [UID](https://en.wikipedia.org/wiki/User_identifier), set a different [Entrypoint](https://docs.docker.com/reference/dockerfile/#entrypoint) than the one specified in the Dockerfile, and other options.
3. A Tool is then turned into an App, which adds the ability to set executable [flag/options](https://unix.stackexchange.com/questions/285575/whats-the-difference-between-a-flag-an-option-and-an-argument), versions, and ultimately is used to execute the software at scale.

> [!NOTE]
> Understandably, this process is not as simplistic as the 3 points above make it sound. In fact, the ingestion process (Docker container > Tool > App) is described in more details here: https://learning.cyverse.org/de/create_apps/#building-an-app-for-your-tool
>
> Additionally, creating a Docker container is not straightforward. Here is some documentation that can help you understand and create your own Docker container:
> - CyVerse's own FOSS content on Containers: [running containers](https://foss.cyverse.org/07_reproducibility_II/) and [building containers](https://foss.cyverse.org/08_reproducibility_III/)
> - [Docker documentation on containerization](https://docs.docker.com/get-started/workshop/02_our_app/)
> - [Docker documentation on running containers](https://docs.docker.com/engine/containers/run/)
> - [Devopscube article on building Docker images](https://devopscube.com/build-docker-image/) (great introductory article showcasing many moving pieces concerning containerization)

### Sharing Data and Apps

Sharing data and Apps with your coworkers can be done using the **Share** function. One can select multiple files (and/or analyses), and click the Share button. This will prompt a small window to open, where the user needs to type the ***exact*** username of the target coworker(s).

Once the usernames are added, users will be able to choose the permission level:

| Permission level | Read | Download/Save | Metadata | Rename | Move | Delete |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Read | X | X | View |			
| Write | X | X | Add/Edit |			
| Own | X | X | Add/Edit | X | X | X |

> [!TIP]
> These are similar to [standard file system permissions](https://en.wikipedia.org/wiki/File-system_permissions), where:
> - Read: allows to read (open/access) the file/folder 
> - Write: allows to modify a file/folder   
> - Own: allows to modify the file/folder and share it with others

#### Teams

<p align="center">
    <img src="https://learning.cyverse.org/assets/de/teams_01.png" width="750">
</p>

As adding each coworker can get tedious, CyVerse allows the creation of <img src="https://learning.cyverse.org/assets/de/menu_items/teamsIcon.svg" width="15"> **Teams**, accessible through the left side menu.

The <img src="https://learning.cyverse.org/assets/de/menu_items/teamsIcon.svg" width="15"> Teams feature allows you to create, organize, and join public or private groups of collaborators, with the goal of enabling a simpler method to share Apps, Tools and data with collaborators.

On the team Teams page, users are able to:

- See all public teams and teams one is part of (top left drop down menu **All Teams**)
- Create a team (top right Team icon)

When creating a Team, one can add as many collaborators as needed and adding a name to the Team.

> [!TIP]
> **REMEMBER TO SAVE BEFORE LEAVING THE TEAMS CREATION PAGE.**

Once finished creating the Team, one can share a file/folder/analysis with the Team instead of each collaborator one by one, thus simplifying the sharing process.

---

## Other Resources: JetStream2

<p align="center">
    <img src="https://docs.jetstream-cloud.org/images/JS2-Logo-Transparent.png" width="350">
</p>

Jetstream2 (JS2) is a user-friendly cloud computing environment for researchers and educators running on [OpenStack](https://www.openstack.org/) and featuring [Exosphere](https://docs.jetstream-cloud.org/ui/exo/exo/) as the primary user interface. It is designed to provide both infrastructure for gateways and other “always on” services as well as giving researchers access to interactive computing and data analysis resources on demand.

JetStream2 is an incredibly powerful resource, that researchers and users can access by submitting a small proposal to [ACCESS](https://access-ci.org/). There are further tiers, however for a smaller research group, the lowest tier is sufficient.

With JetStream2, users can access a variety of powrful CPU nodes, GPU nodes, high memory RAM and large storage nodes. On these nodes, JS2 supports Ubuntu Desktops, JupyterLabs, RStudio, and various Virtual Machine deployments.

Read more on JetStream2 [here](https://docs.jetstream-cloud.org/).

---
---