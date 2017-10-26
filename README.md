K600
====

A GEANT4 simulation of the K600 Spectrometer vault at iThemba LABS, composed of the K600 Spectrometer along with various ancillary detectors 

The primary objective of this simulation is to write the first standard simulation package specifically for the K600 spectrometer that can be utilised by all students/staff that are interested. The compilation of the many different detectors (and other experimental equipment) that are available at iThemba is also of great import. A specific goal to be achieved is to build an accurate simulation of the K600 spectrometer itself by using previous magnet field measurements of the magnets.


I would like to thank the numerous staff at iThemba who put up with my requests to dig through their old documents in the hopes of finding such treasures as detector drawings. Even greater are my thanks to those who have themselves delved into the depths of archives. I hope that this project for iThemba repays your effort.

////////////////////////////////////////////////////////////////////////////////////////////////////

The current implemented detectors include the following:

MMM-type silicon detectors, manufactured by Micron Semiconductor Ltd.
Vertical Drift Chambers (VDC's) which are modelled to the specifications of the newest VDC's at the K600 Spectrometer.
Plastic scintillator detectors known as "paddles" which are used for particle identification and event triggering.
The scattering chamber built for K600, known as BACTAR.
Clover HPGe detectors (Eurogam type), manufactured by Canberra.
The BGO shield complimenting the Clover HPGe detectors (Eurogam type), manufactured by Eurisys.

////////////////////////////////////////////////////////////////////////////////////////////////////

A CAD model import interface called CADMesh (authored primarily by Christopher Poole) is used within this simulation for implementing complex geometries. This needs to be obtained from the following GitHub repository: https://github.com/christopherpoole/CADMesh

Alternatively, one could comment out the relevant CADMesh associated code and use only the hard-coded geometrical objects. It should be noted that an effort has been made to hard-code the geometries in GEANT4 when possible as this has computational advantages.