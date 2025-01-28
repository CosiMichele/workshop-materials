# The moving parts of Functional Open Science

<br>
<br>
<p align="center">
    <img src="https://foss.cyverse.org/assets/foss_title_2024F.jpg" width="600">
</p>
<br>

---
>[!important]
> :clock1: **Schedule**
> - 2:00pm-2:15pm: Introduction to the topic and series (what are we going to learn?) 
> - 2:15pm-2:30pm: What is Open Science?
> - 2:30pm-2:40pm: Version Control
> - 2:40pm-3:00pm: Reproducibility

>[!important]
> :heavy_exclamation_mark: **Requirements**
> - Basic command line knowledge
>- Access to a [Terminal](https://en.wikipedia.org/wiki/Unix_shell)
>    - Unix and Mac users already have access to the Terminal
>    - Windows users can use either [PowerShell](https://en.wikipedia.org/wiki/PowerShell) or the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install)
> - A registered CyVerse account (Register for a CyVerse account)
> - [VSCode](https://code.visualstudio.com/) or any other [Python IDE](https://www.geeksforgeeks.org/top-python-ide/) ( Integrated Development Environment)
> - Some understanding of AI/ML

>[!important]
> :white_check_mark: **Expected Outcomes**
> - Understand what Open Science is
> - Lay the foundations for understanding tools and resources in Open Science
> - Understand the outcomes of the workshop series

<br>

---
---

## Series Overview

**Functional Open Science Skills for AI/ML Applications** is a workshop series that aims to provide attendees with the fundamental knowledge of Open Science and the impact of AI/ML in the field.

As science evolves, scientists are pushed to collaborate in a more transparent fashion (some may even say *pushed* in this direction). Collaborations and transparency are pillars to Open Science, and are becoming requirements with major scientific grants requiring for submissions to have a clear Data/Project Management Plan *before* applying. Additionally, scientific journals such as Nature, PLOS and Science require you to share code and data (if possible) prior to publication.

Thus, it is fundamental for scientists to understand ***what*** Open Science is and ***how*** to best approach today's best practices.

In this workshop series, we aim to educate and introduce the basics of Open Science and AI/ML, reinforcing the learning material through a small functional Open Science friendly AI workflow, applying all that we will be learning through the series. 

---

Today's session focuses on laying the foundations of Open Science by discussing its moving parts. It is important to address that the tools and platforms that we will be discussing are *only a part* of what constitutes Open Science; for a more descriptive dive into Open Science, refere to *the other* [**FOSS**](https://foss.cyverse.org/) -- Foundational Open Science Skills workshop series! FOSS goes deeper into  Open Science, discussing Data Management Plans, Reproducibility, Scalability and more. FOSS will be offerend in the Fall semester, keep your eyes peeled for registrations!

In this session we will cover the following:
- Open Science 
- Version Control 
- Reproducibility
    - Software Environments & Conda
    - Containers

>[!note]
> :question: **Noticing Anything?**
> 
> This session **does not cover the roles of AI/ML**. As these are up and coming technologies, we will be covering them in more depth in the following session.

---
---

## Open Science

<p align="center">
    <img src="https://foss.cyverse.org/assets/open_science_word_cloud.png" width="400">
</p>

Open Science Word Cloud by [*Pownall et al. 2023*](http://dx.doi.org/10.31234/osf.io/vypkb).

The definition of Open Science is not written in stone, making it a topic that may be up to the people to be interpreted. However, there are some attemps at a global definition:

***"Open Science is transparent and accessible knowledge that is shared and developed through collaborative networks"*** - [Vincente-Saez & Martinez-Fuentes 2018](https://linkinghub.elsevier.com/retrieve/pii/S0148296317305441)

or 

***"Open Science is a collaborative and transparent approach to scientific research that emphasizes the accessibility, sharing, and reproducibility of data, methodologies, and findings to foster innovation and inclusivity"*** -ChatGPT

In principle Open Science refers to a set of approaches that would ensure your science to be **transparent**, **collaborative** and **reproducible**.

Obviously, different fields of research have different standardized procedures for experimentation, however oftentimes tools and resources are shared across fields. For example:

- **Need to scale your research?** Scientists can use the [High Performance Computing](https://hpc.arizona.edu/) (HPC) platform or [JetStream2](https://jetstream-cloud.org/index.html) (which we will cover later in the series).
- **Need to share or "save your code"?** [GitHub](https://github.com/) (or [GitLab](https://about.gitlab.com/)) is (are) there for you!
- **Need to make your code reproducible for collaborators (or your future self)?** [Conda](https://docs.anaconda.com/anacondaorg/user-guide/) and [Container technology](https://www.docker.com/) are the current gold standard.

### The Present State of Science

Similar to governmental and other official procedures, the field of science is guided by the outcomes: grants are awarded to groups that have submitted successful proposals, funds are used for research, and the research finds are used to bolster the credibility of the funding agency.

In a sense, ***science is a business***.

Like many businesses, science adapts to the latest trends, and making science Open is latest goal for this cycle.

There are numerous reasons why making science Open is a positive goal:
1. **Accelerates Innovation and Discovery**:
    - Open science allows researchers worldwide to access data, tools, and publications freely. This reduces duplication of effort, fosters collaboration, and speeds up the pace of innovation. Scientists can build on each other’s work more effectively, leading to breakthroughs in less time.
2. **Improves Reproducibility and Credibility**:
    - Sharing data, methodologies, and results openly ensures that studies can be independently verified and reproduced. This strengthens the credibility of scientific findings and reduces instances of irreproducible or questionable research.
3. **Democratizes Knowledge Access**:
    - Open science removes barriers like paywalls, making knowledge accessible to researchers, educators, students, and the general public, regardless of their financial or institutional resources. This promotes equity in learning and research opportunities.
4. **Enhances Public Trust in Science**:
    - When research is transparent and openly shared, it becomes easier for the public to see how conclusions are reached. This fosters trust in scientific processes and can help counter misinformation and skepticism.
5. **Addresses Real-World Problems More Effectively**:
    - Open access to research enables governments, non-profits, and industry to quickly adopt scientific findings to address societal challenges such as climate change, public health crises, and technological advancement. It empowers citizens and organizations to contribute to solutions in meaningful ways.

Furthermore, Open Science (or making science Open) follows a framework that funding agencies are implementing in order for researchers to obtain grants: the Research Life Cycle.

<p align="center">
    <img src="https://foss.cyverse.org/assets/life_cycle.png" width="400">
</p>

The Research Life Cycle from [Open Science Framework](https://osf.io/).

As graduate students, staff and faculty, we are often "forced" to work within a few blocks of the Research Life Cycle. For example, a Data Scientist might be asked to help to collect, store and analyze data whilst a faculty member might be the one driving the research ideas and writing a final report.

We all play a role in the Research Life Cycle, even without us knowing. Often, we find ourselves questioning what we are doing and why we are doing it -- the figure above helps us remind what the big picture is.

>[!note]
> :question: **Last but not least...**
> 
> [... the National Science Foundation (NSF) is pushing to make all publication using taxpayer's money to become open and accessible by 2026.](https://new.nsf.gov/public-access)
> This process requires the use of the framework shown above in the proposals, requiring clear Data Management Plans, detailed list of personnel, collaborations, and material handling.

> **It is of utmost importance that us as scientists understand these requirements to become successful.**

---
---

## The Importance of Version Control

<p align="center">
    <img src="https://swcarpentry.github.io/git-novice/fig/phd101212s.png" width="500">
</p>

We have all been here, taken by the [Software Carpentry Version Control lesson](https://swcarpentry.github.io/git-novice/01-basics.html).

**Version Control** (VC) refers to keeping track of the version of a file, set of files, or a whole project. VC is as much a philosophy as a set of tools; you don't need to master [Git](https://git-scm.com/) to utilize version control (though it is certainly a worthwhile tool for many researchers).

<p align="center">
    <img src="https://content.cdntwrk.com/files/aHViPTg1NDMzJmNtZD1pdGVtZWRpdG9yaW1hZ2UmZmlsZW5hbWU9aXRlbWVkaXRvcmltYWdlXzYzOTkwY2I4OWU5YTUuanBnJnZlcnNpb249MDAwMCZzaWc9OWJjZTA5NDIxNzY4MWFhZjYyNmEwNWNhYmI1YTUzMWQ%253D" width="500">
</p>

The version control path sofware takes before release, from [Cadence](https://resources.pcb.cadence.com/blog/what-is-a-version-control-system).

### What the *Git*?

<p align="center">
    <img src="https://devmountain.com/wp-content/uploads/2022/01/Gitvs_Github-1a-1.jpg" width="500">
</p>

Git vs GitHub, simplified, from [Devmountain](https://devmountain.com/blog/git-vs-github-whats-the-difference/).

**Git** is a command-line program for version control of repositories. It keeps track of changes you make to files in your repository and stores those changes in a .git folder in that repository. These changes happen whenever you make a **commit**. Git stores the history of these commits in a "tree", so you can go back to any previous commit. By keeping track of the **differences** between commits, Git can be much more efficient than storing an entire copy of each version in a document's history.

You could utilize Git completely on its own, on your local computer, and get a lot of benefits. You will have a history of the changes you made to a project, allowing you to go back to any old version of your work. However, where Git really shines is in collaborative work. In order to effectively collaborate with others on a project, you need two basic features: a way to allow people to work in parallel, and a way to host repositories somewhere where everyone can access them. The first feature is **branching**, which is part of Git, and the hosting part can be taken care of by platforms like GitHub, GitLab, or Bitbucket. We will focus on GitHub.

GitHub is a site that can remotely host your Git repositories. By putting your repository onto GitHub, you get a backup of the repository, a way to collaborate with others, and a lot of other features.

- **Git**:
    - First developed in 2005, git is a version control software that allows users to make changes and add versions to their code.
    - Changes and versions are saved locally.
    - Accessible through the Shell.

- **GitHub**:
    - First launched in 2008, its main focus is hosting and sharing code.
    - Uses Git version control software.
    - Changes and versions are saved online (requires an account).
    - Mainly administered through the web (it also has a desktop app).
    - Acquired by Microsoft in 2018.

Here is a great visualization of the Git workflow:

<p align="center">
    <img src="https://www.c-sharpcorner.com/article/git-and-github-version-control-local-and-remote-repository/Images/Git%20And%20Github%20Version%20Control.png" width="500">
</p>

Visualizing the commands through a workflow example
(graphic's correction: ~~marged~~ merged), from [C#Corner](https://www.c-sharpcorner.com/article/git-and-github-version-control-local-and-remote-repository/).

---
---

## Reproducibility

<p align="center">
    <img src="https://foss.cyverse.org/assets/reproducibility-spectrum.png" width="600">
</p>

Source: Peng, *RD Reproducible Research in Computational Science Science* (2011): 1226–1227 via [*Reproducible Science Curriculum*](http://reproducible-science-curriculum.github.io/bosc2015/#/15).

As mentioned throughout today's session, Reproducibility is a key aspect of Open Science. There is no single method to achieve complete reproducibility, as it is a tast that requires diligence, attention and committment. As shown in the graph above, there is a spectrum to Reproducibility. One should *try* to produce a well Reproducible work, but not always it is *necessary* to be on the right side of the spectrum (those are the "unicorn" papers!). Where do you stand on the Reproducibility spectrum?

Thi section is going to cover 3 methods to help you achieve/improve Reproducibility:

- [Software Environments & Conda](#software-environments--conda)
- [Containers](#containers)

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
---
