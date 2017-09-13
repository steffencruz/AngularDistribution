# AngularDistribution

Uses SharcAnalysis and TigressAnalysis libraries to create
angular distributions from ResultsMats histogram file by
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

__A sample input file is shown below__

INFILE: Results_Mats.root  
OUTFILE: Results_1229.root  		  			
CSFILE: AngDist_1229.txt		  					// angular distribution output file    

FRAME:     	Cm  
REACTION:		dp  
NORMALIZATION: 7.93e-4 0.07e-4  

EXC: 				1229

GAM: 		 		415  
GAMPEAK: 		400 425  
GAMBGLO:	 	380 400   
GAMBGHI:	 	425 445  

EXCBINSZ: 	40.0  
GAMBINSZ:   1.0  
THETABINSZ: 4.0  

//bgspec: theta[mid] bglo1 bglo2 bghi1 bghi2 const slope
BGSPEC:			94   			0			0			0			0		 			// REJECT  
BGSPEC:			98   			0			0			0			0		 			// REJECT  
BGSPEC:			102 			0			0			0			0		 			// REJECT  
BGSPEC:			106 			0			0			0			0		 			// REJECT  
