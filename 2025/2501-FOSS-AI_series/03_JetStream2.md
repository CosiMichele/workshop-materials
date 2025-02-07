# Learning to Working in the Cloud: JetStream2 and Reproducibility

<br>
<br>
<p align="center">
    <img src="https://docs.jetstream-cloud.org/images/JS2-Logo-Transparent.png" width="600">
</p>
<br>

---

>[!important]
> :clock1: **Schedule**
> - 2:00pm-2:10pm: Welcome and Introduction to the topic
> - 2:10pm-2:20pm: Getting off your machine (methods for reproducibility)
> - 2:20pm-3:00pm: Introducing JetStream2

>[!important]
> :heavy_exclamation_mark: **Requirements**
> - Basic command line knowledge
>- Access to a [Terminal](https://en.wikipedia.org/wiki/Unix_shell)
>    - Unix and Mac users already have access to the Terminal
>    - Windows users can use either [PowerShell](https://en.wikipedia.org/wiki/PowerShell) or the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install)
> - A registered CyVerse account (Register for a CyVerse account)

>[!important]
> :white_check_mark: **Expected Outcomes**
> - Reproducibility away from your machine (working in the cloud)
> - Knowing where to get compute for research

<br>

---
---

## Getting Off Your Machine (Reproducibility)

Earlier in the workshop series we covered aspects of Reproducibility that would allow not only you but your collaborators to access your software in a much quicker and reliable method. These include creating a [software environment using Conda](https://github.com/ua-datalab/FOSS_AI-ML/wiki/The-moving-parts-of-Functional-Open-Science#software-environments--conda) and [Containers](https://github.com/ua-datalab/FOSS_AI-ML/wiki/The-moving-parts-of-Functional-Open-Science#containers).

Here, we present a summary of these reproducibility methods.

### Software Environments & Conda

Running a pipeline or a workflow constitutes using a number of different programs, tools and software. Often, these software may run only if of a specific version. To ensure that the pipeline runs, one should specify installing the program/software with the right version.

To help with this, we introduce the use of an **Environment Manager** called [**Anaconda** (or **Conda**)](https://docs.conda.io/en/latest/). An Environment Manager allows you to create software installation directories that are isolated from other installations. You can create unique environments and install specific software version to run specific scripts.

Conda is a popular and open source environment manager tool that can be installed on any operating system (Windows, MacOS, Linux).

Users can create environments that have their own set of packages, dependencies, and even their own version of Python.
Another piece of information that's important to remember is that projects can have their own specific requirements, often requiring different versions of the same software. Conda can help you create isolated environments to help softaware not interfering with each other, thus allows for consistent and reproducible results across different systems and setups.

<p align="center">
    <img src="https://miro.medium.com/v2/resize:fit:720/format:webp/0*ElVyaAsDHkIpNgxk.png" width="500">
</p>

Conda, Miniconda, and Anaconda. [Taken from Getting Started with Conda, Medium](https://medium.com/hydroinformatics/getting-started-with-conda-environment-332182d1e937).

---

### Containers

<p align="center">
    <img src="https://foss.cyverse.org/assets/shipping.jpg" width="400">
</p>

Sharing your scientific analysis code with your colleagues is an essential pillar of Open Science that will help push your field forward.

As disscussed, there are technical challenges that may prevent your colleagues from effectively running your code on their computers, issues related to computing environments.

A solution is to package up the code and all of the software and send it to your colleague as a **Container**.

A container is a standard unit of software that packages up code and all its dependencies so the application runs quickly and reliably from one computing environment to another

- Container images are a lightweight, standalone, executable package of software that includes everything needed to run an application: code, runtime, system tools, system libraries and settings
- Each of these elements are specifically versioned and do not change
- The recipient does not need to install the software in the traditional sense

A useful analogy is to think of software containers as shipping containers. It allows us move cargo (software) around the world in standard way. The shipping container can be offloading and executed anywhere, as long the destination has a shipping port (i.e., Docker)

<p align="center">
    <img src="https://cloudblogs.microsoft.com/wp-content/uploads/sites/37/2019/07/Demystifying-containers_image1.png" width="500">
</p>

Difference between Virtual Machines and Containers. Containers are a lot more portable as these do not require an OS to be bundled with the software. Figure source: [Microsoft Cloudblogs](https://cloudblogs.microsoft.com/opensource/2019/07/15/how-to-get-started-containers-docker-kubernetes/).

Containers are similar to virtual machines (VMs), but are smaller and easier to share. A big distinction between Containers and VMs is what is within each environment: VMs require the OS to be present within the image, whilst containers rely on the host OS and the container engine (e.g., Docker Engine).

Software containers, such as those managed by Docker or Singularity, are incredibly useful for reproducible science for several reasons:

- **Environment Consistency**:¶
    - Containers encapsulate the software environment, ensuring that the same versions of software, libraries, and dependencies are used every time, reducing the "it works on my machine" problem.
- **Ease of Sharing**:
  - Containers can be easily shared with other researchers, allowing them to replicate the exact software environment used in a study.
- **Platform Independence**:
  - Containers can run on different operating systems and cloud platforms, allowing for consistency across different hardware and infrastructure.
- **Version Control**:
  - Containers can be versioned, making it easy to keep track of changes in the software environment over time.
- **Scalability**:
  - Containers can be easily scaled and deployed on cloud infrastructure, allowing for reproducible science at scale.
- **Isolation**:
  - Containers isolate the software environment from the host system, reducing the risk of conflicts with other software and ensuring a clean and controlled environment.

<p align="center">
    <img src="https://www.tutorialspoint.com/docker/images/docker_hub_1.jpg" width="600">
</p>

The container's life cycle. Figure source: [Tutorialspoint](https://www.tutorialspoint.com/docker/index.htm).

---

## Working in the Cloud: JetStream2

<p align="center">
    <img src="https://docs.jetstream-cloud.org/images/JS2-Logo-Transparent.png" width="450">
</p>

JetStream2 is a cloud computing environment that runs on the [OpenStack](https://www.openstack.org/) network. It allows to choose various Virtual Machines (VM) sizes with specific requirements, including GPUs (see image below); the only issue I've ran into is the root disk, which maximum size is 60GB. However, you can bypass this problem by attaching a Volume (volumes sizes can go up to TBs of space).

JetStream2's compute is free to access for researchers (trial allocations), however increased use of resources is counted using allocations from ACCESS (formerly known as XSEDE). This means that in order to use JS2 you need an ACCESS account (although ACCESS gives you the possibility to link your institution account, I strongly suggest creating an ACCESS account independent from your institution's account). On top of free access to a small amount of JS2 allocations, ACCESS offers 4 tiers:


| | Description | Allocations awarded | Proposal length | 
|:---:|:---:|:---:|:---:|
| **EXPLORE** | great for resource evaluation, small-medium training events | 400,000 | A paragraph long |
| **DISCOVER** | for modest research needs | 1,500,000 | A page long |
| **ACCELLERATE** | for mid-scale research needs | 3,000,000 | 3 pages long |
| **MAXIMISE** | large scale resource needs | As specified in proposal | Up to 10 pages |

- [content 1](#section-1-subsection-1) 
- [content 2](#section-1-subsection-2)

### Reproducibility

---

### Section 1 Subsection 2

---
---