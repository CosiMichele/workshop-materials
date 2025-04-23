# KEYS Computational Proficiency Assessment 

Welcome to the Computing Skills Assessment. Over the past few weeks, you were asked to complete the following Software Carpentries modules:

- [Shell (also referred to as Command Line or CLI)](https://swcarpentry.github.io/shell-novice/)
- [Python](https://swcarpentry.github.io/python-novice-gapminder/index.html)
- [R](https://swcarpentry.github.io/r-novice-gapminder/)

This is **not a test** of how much you already know — it's a way for us to understand how you work with the basic tools of modern science computing:

You’ll complete **one small script** for each language. Each script should show that you've practiced the **core skills** introduced in the corresponding Carpentries lesson.

We want to see **your thinking and your style** — not a copy-paste answer or something fancy from ChatGPT.

---

## What You’ll Submit

Submit the following **3 scripts**:

| Language | Filename             | Format       |
|----------|----------------------|--------------|
| CLI      | `cli_task.sh`        | Bash script  |
| Python   | `python_task.py`     | Python script|
| R        | `r_task.R`           | R script     |

Each script should:
- Work when run from the terminal.
  - In the case of R, the script can be executed in RStudio.
- Contain **only material found in the Carpentries lesson** for that language.
- Include **comments** explaining your thought process (for example `# The following command appends extracted lines to a new file` or `# remove empty values`).
- Be your own work

---

## Resources

[![](https://raw.githubusercontent.com/CosiMichele/workshop-materials/refs/heads/main/2025/2504-KEYS-assessment/assets/software-carp-overview.png)](https://software-carpentry.org/lessons/index.html)

Each [**Software Carpentry lesson**](https://software-carpentry.org/lessons/index.html) comes with its own data:
- [**Shell/CLI**](https://swcarpentry.github.io/shell-novice/) provides you with [`shell-lesson-data.zip`](https://swcarpentry.github.io/shell-novice/data/shell-lesson-data.zip)
- **Python** is broken down into two parts:
  - [*Programming with Python*](https://swcarpentry.github.io/python-novice-inflammation/) provides:
    - [`python-novice-inflammation-data.zip`](https://swcarpentry.github.io/python-novice-inflammation/data/python-novice-inflammation-data.zip)
    - [`python-novice-inflammation-code.zip`](https://swcarpentry.github.io/python-novice-inflammation/files/code/python-novice-inflammation-code.zip)
  - [*Plotting and programming with Python*](https://swcarpentry.github.io/python-novice-gapminder/) provides:
    - [`python-novice-gapminder-data.zip`](https://swcarpentry.github.io/python-novice-gapminder/files/python-novice-gapminder-data.zip)
- The **R** lesson is also comes in two parts:
  - [*Programming with R*](https://swcarpentry.github.io/r-novice-inflammation/) provides:
    - [`r-novice-inflammation-data.zip`](https://swcarpentry.github.io/r-novice-inflammation/data/r-novice-inflammation-data.zip)
  - [*R for Reproducible Scientific Analysis*](https://swcarpentry.github.io/r-novice-gapminder/) provides:
    - [`gapminder_data.csv`](https://swcarpentry.github.io/r-novice-gapminder/data/gapminder_data.csv)


In case you do are not able to access the files, you can [dowload this zip file (`KEYS-comp-assessment.zip`)](https://github.com/CosiMichele/workshop-materials/raw/refs/heads/main/2025/2504-KEYS-assessment/assets/KEYS-comp-assessment.zip). 

`KEYS-comp-assessment.zip` contains:
- `shell-lesson-data.zip` (needed for the Shell/CLI task)
- `python-novice-gapminder-data.zip` (needed for the Python task)
- `gapminder_data.csv` (needed for the R task)


>[!Note]
> **Folder structure** 
> In order to access the files within, please extract the contents from the zip file first >(how do decompress: [Win](https://support.microsoft.com/en-us/windows/zip-and-unzip-files-8d28fa72-f2f9-712f-67df-f80cf89fd4e5), [Mac](https://support.apple.com/guide/mac-help/zip-and-unzip-files-and-folders-on-mac-mchlp2528/mac), [Ubuntu](stions/86849/how-to-unzip-a-zip-file-from-the-terminal)). Follow this up by extracting the compressed files within as well (`shell-lesson-data.zip` and `python-novice-gapminder-data.zip`) using the same method of decompression.
> 
> - `shell-lesson-data.zip` extracts:
>   ```
>   └── shell-lesson-data/
>       ├── exercise-data/
>       │   ├── alkanes/
>       │   │   ├── cubane.pdb
>       │   │   ├── ethane.pdb
>       │   │   ├── methane.pdb
>       │   │   ├── octane.pdb
>       │   │   ├── pentane.pdb
>       │   │   └── propane.pdb
>       │   ├── animal-counts/
>       │   │   └── animals.csv
>       │   ├── creatures/
>       │   │   ├── basilisk.dat
>       │   │   ├── minotaur.dat
>       │   │   └── unicorn.dat
>       │   ├── numbers.txt
>       │   └── writing/
>       │       ├── LittleWomen.txt
>       │       └── haiku.txt <---------- Use this one!!
>       └── north-pacific-gyre/         
>           ├── NENE01729A.txt
>           ├── NENE01729B.txt
>           ├── NENE01736A.txt
>           ├── NENE01751A.txt
>           ├── NENE01751B.txt
>           ├── NENE01812A.txt
>           ├── NENE01843A.txt
>           ├── NENE01843B.txt
>           ├── NENE01971Z.txt
>           ├── NENE01978A.txt
>           ├── NENE01978B.txt
>           ├── NENE02018B.txt
>           ├── NENE02040A.txt
>           ├── NENE02040B.txt
>           ├── NENE02040Z.txt
>           ├── NENE02043A.txt
>           ├── NENE02043B.txt
>           ├── goodiff.sh
>           └── goostats.sh
>   ```
> - `python-novice-gapminder-data.zip` extracts:
>   ```
>   └── data/
>      ├── gapminder_all.csv
>      ├── gapminder_gdp_africa.csv
>      ├── gapminder_gdp_americas.csv
>      ├── gapminder_gdp_asia.csv
>      ├── gapminder_gdp_europe.csv <---- Use this one!!
>      └── gapminder_gdp_oceania.csv
>   ```

From the extracted files, find:
- `haiku.txt` (`shell-lesson-data/exercise-data/writing/haiku.txt`)
- `gapminder_gdp_europe.csv` (`data/gapminder_gdp_europe.csv`)

For the 3 exercies, you should only need 3 files:

| Comp. Language | Input File | Name of Script |
|:---:|:---:|:---:|
| Shell/CLI | `haiku.txt` | `cli_task.sh` |
| Python | `gapminder_gdp_europe.csv` | `python_task.py` |
| R | `gapminder_data.csv` | `r_task.R` |  


---

## Tasks

### Shell/CLI

Write a script that lists all `.txt` files in a directory and saves the filenames to `txt_files.txt`.

#### Example Script

### Python

 Using only Gapminder data (provided), write a script that prints the average life expectancy for a given continent.

#### Example Script

### R

Read in a CSV, calculate the mean of one column, and plot a simple graph using base R.

#### Example Script

---

## How It Will Be Assessed

Each script will be reviewed using the following criteria:


- Uses only material from lesson
- Script runs without error     
- Output is correct             
- Code is commented/explained   
- Code shows original thinking  

We will **automatically test** (i.e. executed through the Shell (`<script>.sh`, `<script>.py`) and RStudio (`<script>.R`)) your scripts, so please:
- **Name your files correctly**
- **Don’t hardcode file paths or absolute paths**
- **Keep it simple and clean**

---

## A Note on AI

We want to see how **you** think. Submitting code written by AI (ChatGPT, Copilot, Claude, Perplexity, etc.) **goes against the spirit of this assessment**.

This is about **building skill, not getting the answer** — and we’re looking forward to seeing your work!

---

## How to Submit

Please upload your 3 script files in a `.zip` folder to D2L.

