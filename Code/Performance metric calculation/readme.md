# performance_metrics

R code used to obtain dynamic area and dynamic noise of a mean-noise pointset. Briefly, the concave hull of the log10 transformed pointset is converted to a polygon object. Dynamic area is computed as the antilog of the area of the hull. Next, a vertical chord is drawn across the hull at every vertex. Dynamic noise is computed as the antilog of the longest vertical chord. 

## Software used to run code

 R (v3.6.3)
  - ggplot2 (v3.3.0)
  - sp (v1.4-1)
  - alphahull (v2.2)
  - dplyr (v0.8.5)


## Usage
- Run code with a mean-noise pointset as in the test data "mean_noise_data.csv"
- Parameters
   - *df0*: provide path to mean-noise pointset .csv file
   - *alpha*: alpha parameter for concave hull
- Results
  - *F_A*: dynamic area
  - *F_eta*: dynamic noise

## Author
Karl Gerhardt

## License
Copyright (c) 2021, Karl Gerhardt

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
