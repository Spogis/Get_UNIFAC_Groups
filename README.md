# Molecular Fragmentation App for Chemical Engineering

## Overview

The Molecular Fragmentation App is an innovative tool developed to facilitate the a priori calculation of thermophysical properties and the application of predictive thermodynamic models in chemical engineering. It leverages group contribution methods by automating the molecule fragmentation process, which traditionally requires manual effort and extensive databases. This application stands as a significant advancement, offering a rapid and efficient solution for the development and testing of new group contribution methods.

## Features

- **Automated Molecule Fragmentation**: Automates the fragmentation of molecules into chemical groups or other molecular subunits, significantly speeding up the process.
- **Support for UNIFAC Group Contribution Model**: Specifically designed to work with the Universal Quasichemical Functional Group Activity Coefficients (UNIFAC) group contribution model for predicting thermodynamic properties.
- **Database of 20,000+ Molecules**: Tested on a comprehensive database, ensuring wide applicability and reliability.
- **Flexible Heuristic Algorithms**: Implements two heuristic algorithms for fragmentation, accommodating a wide variety of molecular structures.

## Application in Chemical Engineering

The app's capability to automatically fragment molecules is particularly beneficial in chemical engineering for:
- Developing new industrial processes.
- Testing and developing new group contribution methods based on large databases of molecules.
- Applying thermodynamic models or calculating thermophysical properties where predictive methods are necessary.

## Credits

Developed based on the research and models by Simon Müller and described in the paper: 

Müller S. Flexible heuristic algorithm for automatic molecule fragmentation: application to the UNIFAC group contribution model. J Cheminform. 2019 Aug 20;11(1):57. doi: 10.1186/s13321-019-0382-3. PMID: 33430960; PMCID: PMC6701077.

The code utilizes the models developed by the author and is based on the code available at: [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6701077/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6701077/)

## License

MIT License

Copyright (C) 2019, Simon Mueller <simon.mueller@tuhh.de>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## How to Use

1. Install the required dependencies.
2. Run the application and input the molecular structure in SMILES notation.
3. Review the automatically generated fragmentation and the associated thermophysical properties.

For detailed usage instructions and documentation, please refer to the accompanying user guide.
