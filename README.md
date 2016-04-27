# vsGen

Edit Makefile so that the line

CERN_LIBS = -L/usr/lib64/cernlib/2006/lib

makes sense for your system.

Once the Makefile is set up

> make

> ./vsGen -h

The last command will give you instructions on how to use the program.
The output should look like:

SWITCHES:

-h      Print this message

-n<arg> Number of events to generate

-r<arg> User defined random number seed (-r0 = no seeding by user)

-L<arg> Lund output (Lund output text file if -L1, root output file if -L0)

-R<arg> Rxn type:

                0: background (pythia determined. Parameters in pythia.dat)

                1: Omega- K+ K+ K0

                2: Xi- K+ K+

                3: Xi0 K+ K0

                4: Xi-* K+ K+

                5: Xi0* K+ K0

-x<arg> Remove decayed particles if equal to 1

-l<arg> Minimum incident photon energy to generate in GeV: ONLY USED FOR REAL PHOTONS

-u<arg> Maximum incident photon energy to generate in GeV: ONLY USED FOR REAL PHOTONS

-p<arg> Print to screen if equal to 1

-t<arg> Target length in cm

-O<arg> Target center in cm

-v<arg> Vertex resolution in cm

-b<arg> Sigma of beam profile in cm

-o<arg> OutFile name

-e<arg> 0: Real photon. 1: Virtual photon

The above switches overide the default setting.

The default settings are found in the file vsDefaults.conf

The user should modify the vsDefault.conf file to suit their taste

The current default operation is equivalent to the command:

vsGen -n1000 -r0 -R0 -L1 -x0 -l5 -u10.5 -p0 -t7.5 -O0 -v0.1 -b0.01 -e0 -ooutFile.txt
 
