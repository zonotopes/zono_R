# Zonotopes in R

This repository contains a R package for building
the zonotope of a given set of generators and to compute
its volume.

## Getting Started

Clone the git repository recursively:

`git clone --recursive https://github.com/zonotopes/zono_R.git`

## Run demo (from shell)

```sh
Rscript zonotope_demo.R
```

## Run demo (within R)

```
R version 3.4.2 (2017-09-28) -- "Short Summer"
...

> source("zonotope_demo.R")
[1] "N. of generators"
[1] 10
[1] "N. of dimensions (including the output)"
[1] 3
[1] "The volume of the zonotope is"
[1] 20.19476
[1] "Elapsed time (sec):"
elapsed
  0.115
```

## Authors

* **Federico Ponchio** - *Initial work in C++* - [VisLab](http://vcg.isti.cnr.it/~ponchio/)
* **Marco Cococcioni** - *The author of the R porting and owner of this repository* - [Cococcioni @ University Of Pisa](http://www.iet.unipi.it/m.cococcioni/)


## License

This project is granted under simplified BSD license (2-clause).

Copyright (c) 2018, Marco Cococcioni
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the Zonotopes project.

## Acknowledgments

* This work has been funded by professors Marco Grazzi and Giovanni Dosi
* from Scuola Superiore S.Anna, Pisa, Italy
