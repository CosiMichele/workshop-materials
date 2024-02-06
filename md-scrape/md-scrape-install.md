# Installation Journal for MD-Scrape Project

Since the repository for MD-Scrape currently sits on CyVerse, the following software are necessary:

- [GoCommands](https://github.com/cyverse/gocommands)
- [MdRepo](https://github.com/MD-Repo/md-repo-cli)

Users will require CyVerse Accounts in order to use the tools accessing the MD-Scrape Repository.

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
