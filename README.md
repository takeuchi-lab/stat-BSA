# Finding statistically reliable behaviour patterns from hypothesis explosion in data-driven animal behaviour analysis #

This code was used in the experiments in the following papers. For details of the algorithm, please see the following papers.

Finding statistically reliable behaviour patterns from hypothesis explosion in data-driven animal behaviour analysis  
[https://TBA](https://TBA)

Note that this code is created based on S3P-classifier ([https://github.com/takeuchi-lab/S3P-classifier](https://github.com/takeuchi-lab/S3P-classifier)).  

# Verified Environmental

* gcc version 5  
* GNU Make 4.1  

# Setup

## Library
This code uses the following libraries.

   - boost/math/distributions/hypergeometric

Download ([boost libraries](https://www.boost.org/)) and put following location.

   ```
   ./
   ├── include
   │   └── here!
   └── stat-BSA
      └── train
   ```
   
The installation location is specified in the INCLUDE of the Makefile.

## Compile
* make

# Usage
`./train [-options] input_file`

options:  
-    -m : minimum supportSum (default 1)  
    -     If you want to decide on the percentage of the whole, please enter a value less than 1.  
-    -M : minimum length of pattern (default 1)   
-    -L : maximum length of pattern (default 10)  
-    -F : name of reslut file (default output/result.csv)  
-    -p : maximum interval of event (default 0 | -1:none)  
-    -C : whether to do CloSpan (default 0:do not | 1:do)  
-    -R : Repeat resampling in FastWY (default 10000)  
-    -a : significance level (default 0.05)  
-    -s : the Mode of counting supportSum(default 0)  
    -   : 0 is 0 or 1 per record, 1 is the number of pattern  
-    -S : double or single-side test (default 0:double | 1:upper-side | 2:lower-side)  

## Example
`./train -L 100 -F ./output/result.csv ./data/test.txt`

# Demo

## Input data
As shown in the following example, input data should have a label at the beginning and a sequence of symbols behind it.
Note that separate each column with space, and the symbol must a positive integer.  

>1 1 1 1 1 3  
>1 1 1 1 1 3  
>1 1 1 1 1 3  
>1 1 1 1 1 3  
>1 1 1 1 1 2 2  
>1 1 1 1 1 3  
>1 1 1 1 1 3  
>1 1 1 1 1 3  
>1 1 1 1 1 3  
>1 1 1 1 1 2 2  
>1 1 1 1 1 1 3  
>1 1 1 1 1 1 3  
>1 1 1 1 1 1 3  
>1 1 1 1 1 1 3  
>1 1 1 1 1 1 2 2  
>1 1 1 1 1 1 3  
>1 1 1 1 1 1 3  
>1 1 1 1 1 1 3  
>1 1 1 1 1 1 3  
>1 1 1 1 1 1 2 2  
>-1 1 1 1 1 1 2 2 3  
>-1 1 1 1 1 1 2 2 3  
>-1 2 2 3  
>-1 2 2 3  
>-1 2 2 2  
>-1 1 1 1 1 1 2 2 3  
>-1 1 1 1 1 1 2 2 3  
>-1 2 2 3  
>-1 2 2 3  
>-1 2 2 2  
>-1 1 1 1 1 1 2 2 3  
>-1 1 1 1 1 1 2 2 3  
>-1 2 2 3  
>-1 2 2 3  
>-1 2 2 2  
>-1 1 1 1 1 1 2 2 3  
>-1 1 1 1 1 1 2 2 3  
>-1 2 2 3  
>-1 2 2 3  
>-1 2 2 2 

For example, the line  
>1 1 1 1 1 3  

indicates that the label is "1" and the sequence "1 1 1 1 3".

## Output
The output result file shows the adjusted significance level d in the first line.
Each line shows the length, p-value, f-value, total support, support in positive labels, support in negative labels, and which support is more for the extracted significant pattern.

>pattern(d = 0.0138141),length,p,f,supportSum,sup+,sup-,+ or -
>2 2 3,3,1.54171e-07,7.70857e-08,16,0,16,-
>2 3,2,1.54171e-07,7.70857e-08,16,0,16,-
>2,1,1.54171e-07,7.25444e-12,24,4,20,-
>2 2,2,1.54171e-07,7.25444e-12,24,4,20,-
>1 1 1 1 3,5,1.54171e-07,7.70857e-08,16,16,0,+
>1 1 1 3,4,1.54171e-07,7.70857e-08,16,16,0,+
>1 1 3,3,1.54171e-07,7.70857e-08,16,16,0,+
>1 3,2,1.54171e-07,7.70857e-08,16,16,0,+
>1 1 1 1 1 3,6,0.003276,0.001638,8,8,0,+
>...

-----------------------------------------------------------------------
# Licence
GNU General Public License v3.0


