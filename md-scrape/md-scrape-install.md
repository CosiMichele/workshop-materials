# Installation Journal for MD-Scrape Project

Since the repository for MD-Scrape currently sits on CyVerse, the following software are necessary:

- [GoCommands](https://github.com/cyverse/gocommands): allows for movement of CyVerse hosted data 
- [MdRepo](https://github.com/MD-Repo/md-repo-cli): tool for MD data management
- [davfs2](https://wiki.archlinux.org/title/Davfs2): allows for mounting webDAV resources 

Users will require CyVerse accounts in order to use the tools accessing the MD-Scrape Repository.

## Installation of Tools

- Create new folder and add to `PATH`:
  - `sudo mkdir /usr/local/tools && export PATH=$PATH:/usr/local/tools`
- In `/etc/profile.d`, create `md-repo_init.sh`; this will launch the export command everytime the machine is started. The file contains the following lines:
  ```
  #!/bin/bash

  export PATH=$PATH:/usr/local/tools
  ```
- GoCommands:
  - Download GoCommands and decompress: `GOCMD_VER=$(curl -L -s https://raw.githubusercontent.com/cyverse/gocommands/main/VERSION.txt); \ curl -L -s https://github.com/cyverse/gocommands/releases/download/${GOCMD_VER}/gocmd-${GOCMD_VER}-linux-amd64.tar.gz | tar zxvf -`
  - Move to `/usr/local/tools`
- MdRepo:
  - Dowlonad MdRepo and decompress: `CLI_VER=$(curl -L -s https://raw.githubusercontent.com/MD-Repo/md-repo-cli/main/VERSION.txt); \ curl -L -s https://github.com/MD-Repo/md-repo-cli/releases/download/${CLI_VER}/mdrepo-${CLI_VER}-linux-amd64.tar.gz | tar zxvf -`
  - Move to `/usr/local/tools`
- davfs2:
  - installed with `sudo apt install davfs2`
 
## Usage

:exclamation: **Note 1**: As stated above, Users are required to have a CyVerse account. Please go to https://user.cyverse.org/ and create an account. Let me know your username so that you can get access! Access to the folder is restricted to users that have shared their info with me.

:exclamation: **Note 2**: A Google Sheet will be shared with users that need to access the VMs. This Google Sheet will have the username and password that you will be using in order to connect to the VM.

1. Connect to the VM with `<user>@<IP Address>` and paste the password once prompted.
2. Once in the VM, you will need to connect to the CyVerse Data Store. To do so, you will need to type in the following: `sudo mount -o gid=<you>,uid=<you> -t davfs https://data.cyverse.org/dav/iplant/projects/mdrepo-staging/ /opt/md-scrape-data-repo`, where <you> is your username (e.g., `user1`)
3. You should be able to create and remove files using `sudo` commands. This only works with very small files.

:bangbang: **Note 3**: *[With great power comes great responsibility](https://en.wikipedia.org/wiki/With_great_power_comes_great_responsibility)* You have `sudo` priviledges. Before you do something, **THINK ABOUT IT**. :bangbang:
