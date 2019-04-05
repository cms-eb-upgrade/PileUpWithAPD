# PileUpWithAPD

Generate energy deposited in a single EB crystal and APD from N MinBias interactions. APD energies are scaled to "GeV" equivalent as seen by CATIA. Time of APD hits is also generated. It is relative to prompt hits in a crystal and modulo 25 ns

## How to use

There are two steps. First, one need to initialize the channel by choosing its location in eta
````
 EBChannel *ch = new EBChannel(eta);
````

Then, generate energies for each BX with N MinBias interactions
````
 ch->GenerateBX(NPU);
````
Generated energies and times can be accessed as
````
 ch->energySignal()
 ch->energySpike()
 ch->timeSpike()
````

To run an example
````
 root -l -q Example_1EBch_PUandSpikes.C+
````
and inspect histograms and graphs in created
````
 output.root
````
