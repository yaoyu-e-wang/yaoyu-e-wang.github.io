---
title: "Data Wrangling with R Basics - 1"
date: 2020-04-24
draft: false
tags: ["R", "Data Wrangling", "tutorial"]
categories: ["Tutorial"]
author: "Yaoyu E. Wang"

autoCollapseToc: true

---


The R syntax and Data Structre
==============================

R is composed of the following syntax to describe the program: -
*variables* that can be of any of R data types - *operators* operations
that perform a arithematic or numerical action - *functions* a
description of a series of operations to be performed - *comments* the
**\#** indicates comments that are not performed by R

Variables and Data types
------------------------

*Variables* are assignment names for a piece of information, such as
data. The data can take the form of simple *data type* such as integer,
characters, or logics, or can take the form of more complex *data
structure*.

|  Operator | Description                                                    | Example            |
|:---------:|:---------------------------------------------------------------|:-------------------|
|  numeric  | any number including floating                                  | 1, 1.234           |
|  integer  | integer number                                                 | 1,2,3,4…           |
| character | any string of text values declared by quotation **"** or **’** | “a”, “abcd”, ‘123’ |
|  logical  | boolean results                                                | TRUE/FALSE         |

There are other data types that are beyond scope of this course.

We can declare a simple variable and find its data type with the
following:

``` r
# This is a comment that is not read by R 
# Comment is for human to read to explain the code
x=5
x<-10
# Declare a character
y="a"
# Declare a logic
z=TRUE
```

We can also find the data type of a variable by **class** function.

``` r
class(x)
class(y)
class(z)
```

Data Structure
--------------

**Variables** can be stored in different data structure that takes
different types of variable organization.

### Vector

The simplest is the **vector** which contains a series of values with
the same data type.  
We can also declare a vector of numbers into a single variable using
**c()**

For example we can declare the following vectors with different data
types:

``` r
# declare a numeric vector
x=c(1,2,3,4)
# declare a character vector
y=c("a","b","c")
# declare a logic vector, TRUE and FALSE can be short-handed as T and F
z=c(TRUE,FALSE,T,F)

# A data type is assumed when a mixed data types are declared
a=c(1,2,3,"f")
class(a)  # a would be a character vector
```

\[1\] “character”

``` r
b=c(5, F)   # this is a numeric vector  F becomes 0
```

Since **vector** can only be of a single data type, a data type is
assumed when a mixed data types of vector is declared. The default
priority is given to **character**-\>**numeric**/**integer**-\>**logic**

The data type can also be changed into different types:

``` r
class(x)
x=as.character(x)
x
class(x)
x=as.numeric(x)
x
```

We can get the length of the vector as well as access values by calling
for their specific positions using square brackets.

``` r
# length can be determined 
len=length(x)
len
```

\[1\] 4

``` r
x[1]
```

\[1\] 1

``` r
x[2]
```

\[1\] 2

### Matrix

The vector can be tranformed into a matrix by assign a dimension. For
instance, a vector with length of 4 can be transformed into a 2x2 matrix
by **dim** function. A matrix can also be created by **matrix**
function. Again, the data type of the matrix needs to be uniform as
either numeric, character, or logical.

``` r
# change vector into matrix by dim
x=c(1,2,3,4)
dim(x)=c(2,2)
x
```

    ##      [,1] [,2]
    ## [1,]    1    3
    ## [2,]    2    4

``` r
# change vector into matrix by declare matirx
y=c(1,2,3,4)
y=matrix(y, nrow=2, ncol=2)
y
```

    ##      [,1] [,2]
    ## [1,]    1    3
    ## [2,]    2    4

``` r
# declare a character vector into a matrix
a=matrix(c('a','b', 'c', 'd'), nrow=2, ncol=2)
a
```

    ##      [,1] [,2]
    ## [1,] "a"  "c" 
    ## [2,] "b"  "d"

Values within the matrix can be access through defining value
coordinates in very similar way as excel sheet with row and column
numbers.

``` r
x[1,1]
```

    ## [1] 1

``` r
x[2,1]
```

    ## [1] 2

### List

A **list** is similar to vector, but allow for mixing data type of
vectors to be stored in a single variable as elements with assigned
name. A list variable can be declare as a vector of vectors with or
without names. The assigned values still retain its order.

For example:

``` r
# a list variable can be declare as a vector of vectors with or without names
y=list(a=1, 17, b=2:5, c='a')
```

Assigns 4 elements into **y** with the order of 1, 17, 2:5, and ‘a’.
These values can be access by either name with a $ or position within
the list:

``` r
y[[1]]
```

    ## [1] 1

``` r
y[[2]]
```

    ## [1] 17

``` r
y$b
```

    ## [1] 2 3 4 5

We can also get the names from the list variable and rename the
elements:

``` r
names(y)
```

    ## [1] "a" ""  "b" "c"

``` r
names(y)=c("a","b","c","d")
names(y)
```

    ## [1] "a" "b" "c" "d"

``` r
y$b
```

    ## [1] 17

### Data Frame

**Data Frame** is a special kind of list to store rectangular data sets
Think of it as excel sheet where each column is a list.

``` r
df.y=data.frame(y)
df.y
```

    ##   a  b c d
    ## 1 1 17 2 a
    ## 2 1 17 3 a
    ## 3 1 17 4 a
    ## 4 1 17 5 a

The *data.frame* function transform a list into a rectangular data frame
with the number of row is equaled to the longest element. In this case
it is column c, which contains 4 values. All other shorter values are
repeated to fit into the rectangular data frame.

We can use $ to retreive a list within DF, it is equivalent to select a
column in an excel sheet. When we combine list name with a number, it is
equivalent defining a cell in excel with column name and row number like
**‘A1’** points to column A, row 1 cell.

``` r
df.y$a[1]
```

    ## [1] 1

``` r
df.y$c[1]
```

    ## [1] 2

``` r
df.y$c
```

    ## [1] 2 3 4 5

### Factors

**Factor** is a data type unique to R and other statistical language, it
condenses variables into unique categorical variables and store them as
levels. Makes computation of large data with repeated string/character
data type more efficient

``` r
x=sample(letters[1:5], 10, replace = TRUE)
x
```

    ##  [1] "c" "e" "a" "b" "c" "c" "b" "d" "d" "d"

``` r
y=factor(x)
y
```

    ##  [1] c e a b c c b d d d
    ## Levels: a b c d e

``` r
levels(y)
```

    ## [1] "a" "b" "c" "d" "e"

``` r
as.numeric(y)
```

    ##  [1] 3 5 1 2 3 3 2 4 4 4

Fundamental Operators
---------------------

R contains designed arithmetic and logic operators: \#\#\# Arithmetic
\|Operator\|Description\| \|:—:\|:—\| \|+ \|addition \|- \|subtraction
\|\* \|multiplication \|/ \|division \|^ or \*\* \|exponentiation \|x %%
y \|modulus (x mod y) 5%%2 is 1 \|x %/% y \|integer division 5%/%2 is 2

### Logic

|  Operator | Description              |
|:---------:|:-------------------------|
|     \<    | less than                |
|    \<=    | less than or equal to    |
|     \>    | greater than             |
|    \>=    | greater than or equal to |
|     =     | exactly equal to         |
|     !=    | not equal to             |
|     !x    | Not x                    |
|   x \| y  | x OR y                   |
|   x & y   | x AND y                  |
| isTRUE(x) | test if X is TRUE        |

For example we can add two matrices together using **+**

``` r
x=matrix(c(1,2,3,4), nrow=2, ncol=2)
y=matrix(c(5,6,7,8), nrow=2, ncol=2)
x+y
```

    ##      [,1] [,2]
    ## [1,]    6   10
    ## [2,]    8   12

or compare them

``` r
x>y
```

    ##       [,1]  [,2]
    ## [1,] FALSE FALSE
    ## [2,] FALSE FALSE
