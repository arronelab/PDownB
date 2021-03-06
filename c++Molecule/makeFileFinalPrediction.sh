#/bin/sh
rm finalSrc/*.o

g++ -c -O3 -std=gnu++14 -o finalSrc/point.o  finalSrc/point.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o finalSrc/optimizationPoint.o finalSrc/optimizationPoint.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o finalSrc/polyHelix.o finalSrc/polyHelix.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o finalSrc/randomMolGen.o finalSrc/randomMolGen.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o finalSrc/ktlMoleculeRandom.o finalSrc/ktlMoleculeRandom.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o finalSrc/experimentalData.o finalSrc/experimentalData.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o finalSrc/mainPrediction.o finalSrc/mainPrediction.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o finalSrc/hydrationShellRandom.o finalSrc/hydrationShellRandom.cpp -fopenmp
g++ -O3-std=gnu++14 -o predictStructure finalSrc/point.o finalSrc/optimizationPoint.o finalSrc/polyHelix.o finalSrc/randomMolGen.o finalSrc/ktlMoleculeRandom.o finalSrc/experimentalData.o finalSrc/hydrationShellRandom.o finalSrc/mainPrediction.o -fopenmp
