# Data Tapas: Introduction to Markdown 
<br>
<br>
<p align="center">
    <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/4/48/Markdown-mark.svg/1200px-Markdown-mark.svg.png" width="300">
</p>
<br>

---

<br>
<br>

In this Data Tapas workshop we are going to have some fun with **Markdown**!

[Markup languages](https://en.wikipedia.org/wiki/Markup_language) are a system to communicate with your computer how you want a specific word, sentence, table, figure to appear on a website or page. The term "mark up" originates from old manuscripts revisions, where reviewers and editors would "mark up", usually with colour, text to be revised.

<p align="center">
    <img src="https://www.ncbi.nlm.nih.gov/staff/beck/xml/markup/sample.gif" width="600">
    <figcaption style="text-align: center;">
    An "old school" example of marking up. Credits: <a href="your_link_here" style="text-decoration: none; color: blue;">NCBI</a>.
</figcaption>
</p>

Popular markup languages include [LaTeX](https://www.latex-project.org/), [HTML](https://en.wikipedia.org/wiki/HTML), [XML](https://www.w3.org/XML/) and [**Markdown**](https://www.markdownguide.org/basic-syntax/). Notice how many of these are "online" languages, helping to shape a specific website or document.

## Why We Should Care

Markup languages exist to provide a standardized way of annotating and structuring documents or data, making them more readable, interpretable, and accessible to both humans and machines. Additonally, these need to be *lightweight* so that each web page doesn't take a long time to build. 

For example, this website is written in HTML, but the language initially used to tell the computer how to render the text is Markdown. Then the computer is able to read the formatting and present the information in the format we required.  

As part of the effort making science more [Open](https://en.wikipedia.org/wiki/Open_science) and accessible, a lot of the science and collaborations are done online. A lot of code, tutorials, documentation and discussions is now online on platforms such as GitHub, [GitLab](https://about.gitlab.com/), [Stack Overflow](https://stackoverflow.com/), [Read the Docs](https://about.readthedocs.com/?ref=readthedocs.com), [Slack](https://slack.com/), [Reddit](https://old.reddit.com/) and [Jekyll](https://jekyllrb.com/).

Here are some examples of tools/websites (that you may have accessed!) using markdown language:
- [Pandas documentation](https://pandas.pydata.org/docs/): popular data analysis tool; Documentation hosted by read the docs and using an enriched markdown synthax ([ReStructured Text](https://docs.open-mpi.org/en/v5.0.x/developers/rst-for-markdown-expats.html)).
- [PyTorch documentation](https://pytorch.org/docs/stable/index.html): an optimized tensor library for deep learning; Documentation hosted by read the docs and using the same method as above.
- The [DataLab website](https://ua-datalab.github.io/): hosted by [GitHub Pages](https://pages.github.com/) Using [mkdocs-material](https://squidfunk.github.io/mkdocs-material/), an extended markdown language.

Therefore, knowing Markdown will help you with communicating your science, discoveries, and growing your collaborations online!! Because at the end of the day, markdown is communicating with the computer how *you* want your text and websites to look, and therefore allowing potential collaboratiors to understand exactly what you are trying to communicate. **Markdown is a communication tool!**

## Introducing Markdown

Markdown is a [lightweight  markup language (LML)](https://en.wikipedia.org/wiki/Lightweight_markup_language) created to be quickly learned and deployed. 

Released in 2004, Markdown has seen a number of variants rise, such as GitHub Flavoured Markdown (GFM, what you are seeing here), and Markdown Extra (a lot of these features are also seen in GFM). 

Markdown files have the `.md` or `.markdown` extensions.

One of the most useful guides to Markdown is the [markdownguide.org](https://www.markdownguide.org/). There you can find:
- basic synthax guide: https://www.markdownguide.org/basic-syntax/
- Cheat sheets: https://www.markdownguide.org/cheat-sheet/
- Extended synthax (more markdown!): https://www.markdownguide.org/extended-syntax/
- "Hacks" (using HTML for better expression): https://www.markdownguide.org/hacks/
- A list of popular tools that support Markdown editing: https://www.markdownguide.org/tools/

From here on we are going to explore some of the popular Markdown techniques used through the DataLab and beyond! These should render well here on GitHub and other Markdown supporting websites (Stack Overflow, Reddit, Slack, etc...).

### Headings
### Paragraphs
#### Emphasis
##### Bold
##### Italic
##### Bold & Italic
#### Block Quotes
#### Line Breaks
#### Lists
##### Ordered Lists
##### Unordered Lists
#### Images
##### Images with Links
##### Images using HTML
#### Code
##### In-line Code
##### Code Block
#### Links
#### Horizontal Rule
#### Tables

---

## Useful Resources

Markdown can come in various flavours.

- Check out the Markdown guide for basic synthax: https://www.markdownguide.org/basic-syntax/
- Learn about ReStructured Text, useful with Read the Docs: https://docs.open-mpi.org/en/v5.0.x/developers/rst-for-markdown-expats.html
- Write in MkDocs with MkDocs-material for some sleek looking pages: https://squidfunk.github.io/mkdocs-material/reference/

We hope that this workshop has been a useful introduction to markdown for your future communication use!