git clone git@github.com:yguler/DHCALAnalyzer.git

Compile:
g++ `root-config --cflags` -o DHCALAnalyzer TreeAnalysis.cc histDHCAL.cc inputDHCAL.cc anaDHCAL.cc `root-config --glibs`

Run:
./DHCALAnalyzer MuonBeam_Run600009_32GeV 1 0 1


