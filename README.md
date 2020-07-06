# bch
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![GoDoc](https://godoc.org/github.com/zedseven/bch?status.svg)](https://godoc.org/github.com/zedseven/bch)

An implementation of binary Bose-Chaudhuri-Hocquenghem (BCH) codes and error-checking in Go.

## Using the package
To include it in a project, simply use:
```go
import "github.com/zedseven/bch"
```

See [the GoDoc manual](https://godoc.org/github.com/zedseven/bch) for documentation.

## An important note
Please note I did not write the basis of this package - I ported it over from [an excellent example in C on ECCPage](http://www.eccpage.com/bch3.c).
I did, however, clean it up and build the rest of the package around the rudimentary example.
Credit goes to Robert Morelos-Zaragoza for the base code.
