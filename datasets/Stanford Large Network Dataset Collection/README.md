The datasets on this folder where downloaded from the website of the Stanford Network Analysis Project (SNAP), found [here](https://snap.stanford.edu/data/).

More details about the datasets can be found in the table bellow.

| Dataset directory | Description | Nodes | Edges | URL link |
| ----------- | ----------- | ----------- | ----------- | ----------- |
| web-Google | "web-Google" | 875,713 | 5,105,039 | [link](https://snap.stanford.edu/data/web-Google.html) |
| wiki-Talk | "wiki-Talk" | 2,394,385 | 5,021,410 | [link](https://snap.stanford.edu/data/wiki-Talk.html) |

### Adjustments made to the datasets:
The datasets had four (4) lines of meta-data at the beginning of the files and the data were saved in a form that had one edge per line following the pattern "`linkFrom\tlinkTo\n`", like so:
```
linkFrom    linkTo
linkFrom    linkTo
...
```
A program in C was written to discard the meta-data lines and transform the pattern in a new one that has all the out-links of a page in a single row, like so:
```
page_1: linkTo_1 linkTo_2 linkTo_3 ...
page_2: linkTo_1 ...
```

The program is provided in this repository, under the pathname: `/datasets/Stanford Large Network Dataset Collection/graphToAdjacencyList.c`.