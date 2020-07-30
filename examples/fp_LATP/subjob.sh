#!/bin/bash
#=======================================================================
#  Script to be called from fitpot.py and perfom pmd simulation
#  and extract RDF, ADF, and volume data.
#
#  Usage:
#    $ subjob.sh
#=======================================================================

#...copy filed required for pmd calculation
cp ../in.pmd.* ../pmdini ./

#...cd to the directory and clean up
rm -f dump_* out.*

#...Relax pos and cell
# cp in.pmd.relax in.pmd
# pmd 2>%1 > out.pmd.relax
# python ~/src/nap/nappy/pmd/dumps2vol.py --skip=-1 dump_*
# echo "relax done at" `date`

#...NpT MD
cp in.pmd.NpT in.pmd
pmd 2>&1 > out.pmd.NpT
head -n166 out.pmd.NpT
tail -n20 out.pmd.NpT
echo "NpT-MD done at" `date`
#...extract rdf, adf, vol and rename files
python ~/src/nap/nappy/rdf.py -d 0.05 -r 5.0 --gsmear=2 --skip=80 --specorder=Li,Al,Ti,P,O --pairs Li-O,Al-O,Ti-O,P-O,O-O --out4fp -o data.pmd.rdf pmd_[0-9]* 2>&1 
python ~/src/nap/nappy/adf.py --gsmear=2 --triplets=Al-O-O,Ti-O-O,P-O-O --skip=80 --out4fp -o data.pmd.adf pmd_[0-9]* 2>&1 
python ~/src/nap/nappy/vol_lat.py --skip=80 --out4fp pmd_[0-9]* 2>&1
#tail -n1 out.erg | awk '{print $7}' > data.pmd.vol
echo "post-processing done at" `date`
