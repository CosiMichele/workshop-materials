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

Headinds are created by adding `#` in front of a word or a phrase. The more the `#` the "smaller" the level heading: 1 `#` represents a level one heading used for the name of a chapter, 2 `##` is a level 2 heading used, for example, as heading of a paragraph, all the way until level heading 6. 

|Markdown|Output|
|---|---|
|# Heading level 1 | <h1>Heading level 1</h1> |
|## Heading level 2 | <h2>Heading level 2</h2> |
|### Heading level 3 | <h3>Heading level 3</h3> |
|#### Heading level 4 | <h4>Heading level 4</h4> |
|##### Heading level 5 | <h5>Heading level 5</h5> |
|###### Heading level 6 | <h6>Heading level 6</h6> |

### Paragraphs and Line Breaks

Paragraphs can be populated as normal text, however you will require at least 2 line breaks (or character returns) in order to add separation between 2 sentences.

Wrong Markdown:
```
This is my first paragraph. 
This is my second paragraph. 
There's no separation as I have only 1 character return.
```

Output: 

This is my first paragraph. 
This is my second paragraph. 
There's no separation as I have only 1 character return.

Right Markdown:
```
This is my first paragraph. 

This is my second paragraph. 

There is separation as I have put 2 character returns between the two sentences.
```

Output:

This is my first paragraph. 

This is my second paragraph. 

There is separation as I have put 2 character returns between the two sentences.

#### Emphasis

Markdown allows you to to emphasize sentences or words in **bold**, *italic* or ***both***!

||Markdown|Output|
|---|---|---|
|Bold| `**Use asterisks twice around the word or sentence to bold it!**` | **Use asterisks twice around the word or sentence to bold it!** |
|Italic| `*Use asterisk once round the word or sentence to italicize it!*` | *Use asterisk once round the word or sentence to italicize it!* |
|Bold and Italic | `***Use 3 asterisks around the word or sentence to have both bold and italic!***` | ***Use 3 asterisks around the word or sentence to have both bold and italic!*** |

#### Block Quotes

Block quotes can be used in order to emphasize a note or a warning. These are created using the `>` character.

Markdown:

`> Note: 'tis a block quote!`

Output:

> Note: 'tis a block quote!

There can be multiple sentences in block quotes. Simply remember to add a `>` even for emtpy lines.

Markdown:

```
> This is the first sentence in the block quote...
>
> ... and here's the second!
```

Output:

> This is the first sentence in the block quote...
>
> ... and here's the second!

You can also create nested blocked quotes by adding more `>` in sequence.

Markdown:

```
> This is the heading sentence.
>
>> Here's a nested sentence...
>>> And here's a veeeery (x3) nested sentence!
```

Output:

> This is the heading sentence.
>
>> Here's a nested sentence...
>>> And here's a veeery (x3) nested sentence!

#### Lists

Lists can either be ordered or unordered.

##### Ordered lists

Ordered lists can be called with simply typing `1.` in the first sentence and the list will continue from there automatically. Even if the second number is NOT 2, the list will resume from 2. If the number starts from a different number, the list will start that specified number.

Inserts however, need to start from `1.`

Markdown:
```
1. Sentence 1
1. Sentence 2 
    1. Nested insert 1
    1. Nested insert 2
1. Sentence 3

Break

4. Sentence 1
4. Sentence 2
    1. Nested insert 1
    1. Nested insert 2
4. Sentence 3
```

Output:

1. Sentence 1
1. Sentence 2 
    1. Nested insert 1
    5. Nested insert 2
1. Sentence 3

Break

4. Sentence 1
4. Sentence 2
    1. Nested insert 1
    5. Nested insert 2
4. Sentence 3
    
##### Unordered Lists

Unordered lists are called with `-`, `*`, `+`.

Markdown:
```
- Sentence 1
- Sentence 2 
    - Nested insert 1
    - Nested insert 2
- Sentence 3

Break

+ Sentence 1
+ Sentence 2 
    + Nested insert 1
    + Nested insert 2
+ Sentence 3

Break

* Sentence 1
* Sentence 2 
    * Nested insert 1
    * Nested insert 2
* Sentence 3
```

Output:

- Sentence 1
- Sentence 2 
    - Nested insert 1
    - Nested insert 2
- Sentence 3

Break

+ Sentence 1
+ Sentence 2 
    + Nested insert 1
    + Nested insert 2
+ Sentence 3

Break

* Sentence 1
* Sentence 2 
    * Nested insert 1
    * Nested insert 2
* Sentence 3

#### Images

Images are insterted using the following convenctions:

Markdown:
```
![image name](URL-to-image)
[![image name](URL-to-image)](URL destination)

such as 

![wilma&wilbur](https://alumni.arizona.edu/sites/default/files/styles/az_natural/public/2022-06/wilbur%20and%20wilma.jpg?itok=BtBNnOAR)
[![wilma&wilbur](https://alumni.arizona.edu/sites/default/files/styles/az_natural/public/2022-06/wilbur%20and%20wilma.jpg?itok=BtBNnOAR)](https://alumni.arizona.edu/history-traditions/wilbur-and-wilma)
```

Output:

The first image will just show the image.
![wilma&wilbur](https://alumni.arizona.edu/sites/default/files/styles/az_natural/public/2022-06/wilbur%20and%20wilma.jpg?itok=BtBNnOAR)

The second picture will redirect to a selected webpage.
[![wilma&wilbur](https://alumni.arizona.edu/sites/default/files/styles/az_natural/public/2022-06/wilbur%20and%20wilma.jpg?itok=BtBNnOAR)](https://alumni.arizona.edu/history-traditions/wilbur-and-wilma)

#### Code

Code can be shown in 2 ways, either in-line or in a code block.

##### In-line Code

When addressing in-line code, all you need is adding ` (above the ~ sign) next to your word or sentence.

Markdown:
```
`This sentence is going to render as code.`
```

Output:

`This sentence is going to render as code.`

##### Code Block

The Code Block synthax is what we've been using throughout the page in order to show you what the raw Markdown is! In order to create code blocks, one has to add 3 ticks at the beginning and at the end of a code chunk.

Markdown:
````
```
This is the first sentence in a code block

code line 1:
    code line 2
        code line 3

This is going to be the last sentence in the code block.
```
````

Output:

```
This is the first sentence in a code block

code line 1:
    code line 2
        code line 3

This is going to be the last sentence in the code block.
```

#### Links

We have seen how to link things earlier during the images section and we can apply the same synthax to words and sentences. One has to follow the next synthax: `[sentence in square brakets](followed by the URL destination).`

Mardown:

```
[This sentence is going to be linking to the Data Science Institute webpage.](https://datascience.arizona.edu/)
```

Output:

[This sentence is going to be linking to the Data Science Institute webpage.](https://datascience.arizona.edu/)


#### Horizontal Rule

Adding orizontal separators in a page is simple: add 3 dashes (`-`) to create a visual separation!

Markdown:
`---`

Output:

---

#### Tables

Tables are created by using pipes (`|`) and dashes (`-`).

Markdown:

```
|Entry A|Entry B| Data |Results|
|---|---|---|---|
|entry 1A|entry 1B|data 1|result 1|
|entry 2A|entry 2B|data 2|result 2|
|entry 3A|entry 3B|data 3|result 3|
|entry 4A|entry 4B|data 4|result 4|
```
Output:

|Entry A|Entry B| Data |Results|
|---|---|---|---|
|entry 1A|entry 1B|data 1|result 1|
|entry 2A|entry 2B|data 2|result 2|
|entry 3A|entry 3B|data 3|result 3|
|entry 4A|entry 4B|data 4|result 4|

Notice how the first line is always the table header and the second always a separator. You can add handendness in the second line to set how that specific column is going to align using `:` on the left, right, or both for centered. 

You can also leave columns (or rows) blank to create space.

Markdown:

```
|Left handed column|Centered Column| Empty column | Right handed column |
|:---|:---:|---|---:|
|entry 1A|entry 1B||data 1|
|entry 2A|entry 2B||data 2|
|empty row||||
|entry 3A|entry 3B||data 3|
|entry 4A|entry 4B||data 4|
```
Output:

|Left handed column|Centered Column| Empty column | Right handed column |
|:---|:---:|---|---:|
|entry 1A|entry 1B||data 1|
|entry 2A|entry 2B||data 2|
|empty row||||
|entry 3A|entry 3B||data 3|
|entry 4A|entry 4B||data 4|


---

## Useful Resources

Markdown can come in various flavours.

- Check out the Markdown guide for basic synthax: https://www.markdownguide.org/basic-syntax/
- Learn about ReStructured Text, useful with Read the Docs: https://docs.open-mpi.org/en/v5.0.x/developers/rst-for-markdown-expats.html
- Write in MkDocs with MkDocs-material for some sleek looking pages: https://squidfunk.github.io/mkdocs-material/reference/

We hope that this workshop has been a useful introduction to markdown for your future communication use!