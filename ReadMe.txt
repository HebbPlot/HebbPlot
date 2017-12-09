Requirements:

Perl (https://www.perl.org)
Bedtools (http://bedtools.readthedocs.io/en/latest/)
Matlab (https://www.mathworks.com)
Statistics and Machine Learning Toolbox (https://www.mathworks.com/products/statistics.html)
PDF Viewer (https://get.adobe.com/reader/)
UNIX/Linux/Mac 

To install, run the following command:
./install.pl

To learn about the paramaters of HebbPlot, run the following command:
./hebbPlot.pl

Note about sort and merge commands: if your epigenomes aren't sorted/merged, you can tell hebbPlot
to sort and merge them for you. After doing it once, you can reuse those same sorted/merged epigenomes.
Simply move the overlap/ directory to and specify this directory the next time you run hebbPlot again.


If you want to produce results similar to Table 1, as seen in the paper,
use comparePosNegProm.m in the Matlab directory.

If you want to produce results similar to Table 2 and 3, as seen in the paper,
use averageDotSim.m in the Matlab directory.