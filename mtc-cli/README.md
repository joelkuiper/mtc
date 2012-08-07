GeMTC CLI - Mixed Treatment Comparison model generation
=======================================================

[GeMTC][1] CLI is a program for Mixed Treatment Comparison (MTC) model
generation. MTC is also known as Network Meta-Analysis.

The software is capable of generating JAGS and BUGS models for both
dichotomous and continuous data structures. The data is provided in an
XML file and processed by the library to generate the correct
parametrization, whether a consistency, an inconsistency or a node-split
model.

Several flavours of GeMTC are available:

  * A graphical user interface (GUI) application with built-in data
    management and model generation capabilities. The generated models
    can be run using JAGS or BUGS.

  * A command line interface (CLI) program to generate JAGS/BUGS from
	data files provided on the command line. The generated models can be
    run using JAGS or BUGS.

  * An R package that combines the GeMTC library and JAGS to provide
    full analysis capabilities (EXPERIMENTAL).

  * A Java library  that can also run the MTC models using the YADAS
    pure Java MCMC package. This option is used by [ADDIS][2] to perform
    network meta-analysis.

Several example XML files are provided with the GeMTC distribution to
demonstrate the capabilities of the library. All examples were taken
from published network meta-analyses, and the references to the original
papers are in the respective XML files.

Requirements
------------

GeMTC is built in [Java][3], and requires Java 6 (JRE 1.6) or newer. To
run the generated models, you need [JAGS][4] 3.0 or newer, or
[OpenBUGS][5], or [WinBUGS][6]. You may want to use [R][7] to perform
post-processing and analysis.

Using GeMTC CLI with JAGS
-------------------------

For example, say we have a data file 'cardiovasc.xml'. To generate and
run the model, do the following (on a command line; in windows, use
start -> run -> type 'cmd' -> ok). Place the 'mtc.R' file provided in
with MTC in the same directory as your data.

 $ cd directory_with_data_file_and_mtc

 $ java -jar mtc-0.10.jar --type=consistency cardiovasc.xml cv.cons

    ... some output follows ...   

 $ jags cv.cons.script

    ... running model ...

 $ R

 > source('mtc.R') 
 > data <- read.mtc('cv.cons')
 > summary(data)

    ... summary of simulation results ...

To run an inconsistency model, replace --type=consistency with
--type=inconsistency. The following are important ways to assess
convergence:

 > plot(data)
 > gelman.diag(data)
 > gelman.plot(data)

For consistency and inconsistency models, a utility method is provided
to calculate values for derived (functional) parameters. This is used
like this:

 > source('mtc.R') 
 > source('cv.cons.analysis.R')
 > data <- append.derived(read.mtc('cv.cons'), deriv)
 > summary(data)

    ... summary of simulation results ...

Node split models can be generated with --type=nodesplit, and utility
methods are also provided in "mtc.R":

 > source('mtc.R') 
 > data <- read.mtc('cv.splt.A.B')
 > summary(data)

    ... summary of simulation results ...

 > nodeSplitP(data, 'A', 'B')

    ... P-value for A.B.ind == A.B.dir

To evaluate the P-values for all generated node-split models, use:

 > nodeSplitSummary('cv')

Versions
--------

0.12.3: Bugfixes.

0.12.2: Bugfixes and minor improvements.

0.12.1: Sort studies and treatments alphabetically, fix flaw in prior
generation, and other bugfixes. Improve integration with ADDIS and allow
the simulation to be extended after assessment of convergence.

0.12: Complete rewrite in Java, redesigned node-splitting model
generation.

0.10.1: Bugfixes, generate R code for JAGS models.

0.10: Addition of GUI module, improved consistency model generation,
support for BUGS. Renamed to GeMTC CLI and GeMTC GUI.

0.8: Starting value generation for JAGS models, JAGS node-split models,
provide R code for analysis.

0.6: Node-splitting models, assessment of convergence, starting value
generation, allow access to samples (all for YADAS only, JAGS models
have not advanced). Modularized into 4 libraries.

0.4: Improved parametrization algorithm. Allow output of Network as XML.
Check connectedness of network. Some issues fixed.

0.2: Initial release. JAGS and YADAS model generation. Java API for
running MTC models.

Building from source
--------------------

Source code can be obtained from [the GeMTC site][1]. In order to build
GeMTC you need:

 - [Java >= 1.5][3]
 - [Maven 2][8]
 - Other dependencies are downloaded automatically by Maven

To build, use "mvn package". Due to downloading dependencies, the first
time this may take a while.

License
-------

GeMTC is open source, and licensed under GPLv3. See LICENSE.txt for more
information.

Contact
-------

Contact the author, Gert van Valkenhoef, g.h.m.van.valkenhoef@rug.nl,
for more information.


[1]: http://drugis.org/gemtc
[2]: http://drugis.org/addis
[3]: http://www.java.com/getjava/
[4]: http://sourceforge.net/projects/mcmc-jags/
[5]: http://www.openbugs.info/
[6]: http://www.mrc-bsu.cam.ac.uk/bugs/
[7]: http://r-project.org/
[8]: http://maven.apache.org/
