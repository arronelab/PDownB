#!/bin/bash
ScatterFile=newFitData/$1/$1Saxs.dat
SequenceFile=newFitData/$1/fingerPrint.dat
initialCoordsFile=newFitData/$1/coordinates.dat
pairedPredictions=none
fixedsections=newFitData/$1/varyingSectionSecondary.dat
crystalSymmetry=none
withinMonomerHydroCover=none
betweenMonomerHydroCover=none
kmin=0.022;
kmax=0.25;
NoSteps=10

mkdir newFitData/$1/$2


for i in {1..2}
do
    nChanges $ScatterFile $SequenceFile $initialCoordsFile $pairedPredictions $fixedsections $crystalSymmetry $withinMonomerHydroCover $betweenMonomerHydroCover $kmin $kmax $NoSteps newFitData/$1/$2/mol
done
 
