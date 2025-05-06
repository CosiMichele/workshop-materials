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

| Language | Filename | Format |
|:---:|:---:|:---:|
| :shell: Shell/CLI | **`cli_task.sh`** | Bash script |
| :snake: Python | **`python_task.py`** | Python script |
| :registered: R | **`r_task.R`** | R script |

Each script should:
- Work when run from the terminal.
  - In the case of R, the script can be executed in RStudio.
- Contain **only material found in the Carpentries lesson** for that language.
- Include **comments** explaining your thought process (for example `# The following command appends extracted lines to a new file` or `# remove empty values`).
- Be your own work

---

## :globe_with_meridians: Resources

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
> In order to access the files within, please extract the contents from the zip file first >(how do decompress: [Win](https://support.microsoft.com/en-us/windows/zip-and-unzip-files-8d28fa72-f2f9-712f-67df-f80cf89fd4e5), [Mac](https://support.apple.com/guide/mac-help/zip-and-unzip-files-and-folders-on-mac-mchlp2528/mac), [Ubuntu](https://askubuntu.com/questions/86849/how-to-unzip-a-zip-file-from-the-terminal)). Follow this up by extracting the compressed files within as well (`shell-lesson-data.zip` and `python-novice-gapminder-data.zip`) using the same method of decompression.
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
- `haiku.txt` (found in `shell-lesson-data/exercise-data/writing/haiku.txt`)
- `gapminder_gdp_europe.csv` (found in `data/gapminder_gdp_europe.csv`)

For the 3 exercies, you should only need 3 files:

| Comp. Language | Input File | Name of Script to Submit |
|:---:|:---:|:---:|
| :shell: Shell/CLI | **`haiku.txt`** | `cli_task.sh` |
| :snake: Python | **`gapminder_gdp_europe.csv`** | `python_task.py` |
| :registered: R | **`gapminder_data.csv`** | `r_task.R` |  


---

## :notebook: Tasks

### :shell: Shell/CLI

For the Shell exercise, use the `haiku.txt` file obtainable from [`shell-lesson-data.zip`](https://swcarpentry.github.io/shell-novice/data/shell-lesson-data.zip), within the [Shell/CLI](https://swcarpentry.github.io/shell-novice/) lesson.

<u>

**Write a script that:** 

1. **Finds and navigates to the `haiku.txt` file**
2. **Print its location to the terminal**
3. **Count the words of the file**
4. **Append data at the end of the file and print the final product**

</u>

This task is completely doable through material found in the [Shell/CLI](https://swcarpentry.github.io/shell-novice/) carpentry lesson, with ~10 lines.

#### :shell: Example Task and Script

In this example task, we extract lines with a keyword ("not"), count matching words, append summary, and print the file. It covers navigation (slightly), finding files and searching within the file and appending. You may find yourself using some of these commands!

```bash
#!/bin/bash

# Find haiku.txt starting from current directory
file_path=$(find . -name "haiku.txt" | head -n 1)

# Navigate to its directory
cd "$(dirname "$file_path")"

# Print absolute path
echo "haiku.txt found at: $(pwd)/haiku.txt"

# Define keyword to search
keyword="not"

# Count how many times the keyword appears
keyword_count=$(grep -o "$keyword" haiku.txt | wc -l)

# Extract and append keyword context
echo "" >> haiku.txt
echo "---- Keyword Summary ----" >> haiku.txt
echo "Keyword '$keyword' appears $keyword_count times" >> haiku.txt
echo "Summary generated on: $(date)" >> haiku.txt
echo "--------------------------" >> haiku.txt

# Show updated file
cat haiku.txt
```

---

### :snake: Python

 Using only Europe Gapminder data (`data/gapminder_gdp_europe.csv`), <u>**write a script that lists all countries where GDP in 1952 was less than 5000 and are above 5000 in 2007. Plot the GDP per capita over time.**</u>

This task should take roughly 20 lines of code, with ~5 are needed for the plotting and 2 for importing.

**NOTE:** make sure that this code is executable from the parent folder of `data/` (e.g., `Documents/data/gapminder_gdp_europe.csv`).

#### :snake: Python Example Task and Script

In this example, we are asked to plot the GDP per capita over time for selected countries. We use 15 lines of code, 2 of which are loading libraries, 5 of which are plotting-related. All the commands used are present in the [Plotting and programming with Python](https://swcarpentry.github.io/python-novice-gapminder/) Software Carpentry page.

Feel free to execute the code yourself!

```python

# Load libraries
import pandas as pd
import matplotlib.pyplot as plt

# Load the dataset
data = pd.read_csv('data/gapminder_gdp_europe.csv')

# Set of countries to plot
countries_to_plot = ['Germany', 'France', 'Italy', 'Spain']

# Loop through each country and plot its data
for country in countries_to_plot:
    country_data = data[data['country'] == country]
    # Extract year labels and values (skip the 'country' column)
    years = country_data.columns[1:]
    values = country_data.iloc[0, 1:]
    plt.plot(years, values, label=country)

# Add labels and title
plt.xlabel('Year')
plt.ylabel('GDP per Capita')
plt.title('GDP Per Capita Over Time')
plt.legend()
plt.xticks(rotation=45)

# Show the plot
plt.show()
```

---

### :registered: R

 Using the Gapminder data provided by [R for Reproducible Scientific Analysis](https://swcarpentry.github.io/r-novice-gapminder/) ([`gapminder_data.csv`](https://swcarpentry.github.io/r-novice-gapminder/data/gapminder_data.csv)), <u>**write a script that plots a line graph of the top 10 countries with lowest life expectanct in 1952, and their respective life expectancy progress over time**</u>

This task should take ~15 lines of code. All of the code required is accessible through the [R for Reproducible Scientific Analysis](https://swcarpentry.github.io/r-novice-gapminder/) Software Carpentry lesson.

Ensure that the code can be viewed and executed in RStudio using `gapminder_data.csv`.

#### :registered: R Example Task and Script

In this example, we are creating 2 boxplots comparing the GDP per capita of the African countries vs European countries in the year 2007. Hint: some of the functions used here will be useful in the assigned task.

```r
# Load tidyverse
library(tidyverse)
setwd("~/")

# Read data
gapminder <- read.csv("gapminder_data.csv")

# Filter for data in 2007 for two selected continents
gapminder_2007 <- gapminder %>%
  filter(year == 2007, continent %in% c("Africa", "Europe"))

# Plot GDP per capita comparison using boxplot
ggplot(data = gapminder_2007, aes(x = continent, y = gdpPercap)) +
  geom_boxplot() +
  labs(title = "GDP per Capita in 2007: Africa vs. Europe",
       x = "Continent", y = "GDP per Capita")

```

---

## :white_check_mark: How It Will Be Assessed

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

## :heavy_exclamation_mark: A Note on AI

We want to see how **you** think. Submitting code written by AI (ChatGPT, Copilot, Claude, Perplexity, etc.) **goes against the spirit of this assessment**.

This is about **building skill, not getting the answer** — and we’re looking forward to seeing your work!

---

## How to Submit

Please upload your 3 script files individually or in a compressed `.zip` folder to D2L.

