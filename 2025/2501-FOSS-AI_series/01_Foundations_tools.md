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
    - Workflows

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

- [Software Environments & Conda](#software-environments--conda)
- [Containers](#containers)
- [Workflows](#workflows)

### Software Environments & Conda

---

### Containers

<p align="center">
    <img src="https://foss.cyverse.org/assets/shipping.jpg" width="400">
</p>

<p align="center">
    <img src="https://cloudblogs.microsoft.com/wp-content/uploads/sites/37/2019/07/Demystifying-containers_image1.png" width="500">
</p>

Difference between Virtual Machines and Containers. Containers are a lot more portable as these do not require an OS to be bundled with the software. Figure source: [Microsoft Cloudblogs](https://cloudblogs.microsoft.com/opensource/2019/07/15/how-to-get-started-containers-docker-kubernetes/).

<p align="center">
    <img src="https://www.tutorialspoint.com/docker/images/docker_hub_1.jpg" width="600">
</p>

The container's life cycle. Figure source: [Tutorialspoint](https://www.tutorialspoint.com/docker/index.htm).

---

### Workflows


<p align="center">
    <img src="https://raw.githubusercontent.com/nf-core/rnaseq/3.15.1/docs/images/nf-core-rnaseq_metro_map_grey_animated.svg" width="600">
</p>

Image source: [nf-core/rnaseq](https://nf-co.re/rnaseq/3.15.1/).

---
---

## Bringing It All Together