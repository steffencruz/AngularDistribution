# AngularDistribution

Uses SharcAnalysis and TigressAnalysis libraries to create
angular distributions from RedwoodMats histogram file by
making projections of an ExcVsTheta matrix and extracting
the counts within a specified excitation energy window

___Input Description___

A standardized input file is required to create an angular 
distribution. This describes all energy and angular gates 
to be used, including background shapes. The input can also
be specified and included using the CNTSPEC option, which is
useful for elastic channels when the peaks overlap.

The full list of accepted input specifications can be found 
in the structure entitled ' SteffenOptions '

An efficiency correction is automatically carried out if available
as is a branching ratio correction, see TTigressAnalysis docs.
