# Graphite

GRAPHITE
GRAph used tool to correct, PHase, and ImpuTE VCF


## Usage

graphite  -i VCF [--ref ref VCF] -p ped [-m map] [-t num_threads] [-f family indices] [-c chrom indices] [--lower-progs lower num progenies] [--large-only] [--not-impute-isolated [--out-isolated]] [--not-correct-isolated] -o out.
family indices: (index|first:last)[,(index|first:last)[,..]]
chrom indices: same as family indices.

*VCF*			  : input VCF

*ped*			  : pedigree file
					A text file with 4 columns and without header.
					Its delimiter is a space.
					Its columns are family name, sample name, a parent name, and other parent name.
					Parent order is reflected in the VCF's Genotype.

*out*			  : output VCF

### optional

*ref VCF*		  : reference VCF

*map*			  : genetic map file
					A CSV file with 3 columns and without header.
					Its deliimiter is a comma.
					Its columns are scaffold name, cM, and Mbp.
					If not specified, 1cM = 1Mbp.

*chrom indices*   : output scaffolds' indices

*lower num progenies* : the lower limit on the number of progeny that can be regarded as large families.

--large-only	  : output large families only

--not-impute-isolated : not impute isolated samples

--out-isolated	  : output not imputed isolated samples

--not-correct-isolated : impute and phase, but not correct isolated samples

#### 

## License
MIT License

Copyright (c) 2023 Minoru Inamori

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
