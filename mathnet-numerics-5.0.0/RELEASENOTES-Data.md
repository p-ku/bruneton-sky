### 4.1.0 - 2020-12-30
* Support for Math.NET Numerics v4.12
* Resources inlined *~Philipp Huber*

### 4.0.0 - 2018-02-14
* Support for Math.NET Numerics v4
* In addition to .NET 4.0 and newer now also targets .NET Standard 1.3 and 2.0.

### 4.0.0-beta01 - 2018-02-04
* Support for Math.NET Numerics v4
* In addition to .NET 4.0 and newer now also targets .NET Standard 1.3 and 2.0.

### 3.2.1 - 2017-04-29
* BUG: to prevent corrupt files, writers now overwrite files if they exists already

### 3.2.0 - 2016-04-11
* Text: reduce memory usage of parsing *~liuzh*

### 3.1.1 - 2015-07-13
* MatrixMarket: fix MatrixMarketWriter.WriteVector

### 3.1.0 - 2015-01-11
* Text: support for missing values in delimited reader *~Marcus Cuda*

### 3.0.0 - 2014-07-23
* Requires at least Math.NET Numerics v3.0.0
* Text: static-only DelimitedReader/Writer api
* Matlab: static/stateless writer api design, consistent with reader
* Matlab: much more efficient reading/writing of sparse data *~Christian Woltering*
* Matlab: intermediate MatlabMatrix objects to support mixed-type matrices

### 3.0.0-beta02 - 2014-06-15
* Matlab: MATLAB Reader now static and stateless; easier to use
* Text: MatrixMarket: Compression enum instead of boolean flag
* Require Math.NET Numerics v3.0.0-beta03 (fixes build-version issue)

### 3.0.0-beta01 - 2014-04-23
* First v3 beta release
