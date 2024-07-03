
#python coverage.py -s ../data/GW190425/LALInference.fits.gz -f ../data/GW190425/ZTF_fields.dat -l ../data/GW190425/transients.dat -o ../output/GW190425/ZTF_LALInf.pdf -t ZTF

#python coverage.py -s ../data/GW190425/bayestar.fits.gz -f ../data/GW190425/ZTF_fields_bayestar.dat -l ../data/GW190425/transients.dat -o ../output/GW190425/ZTF_bayestar.pdf -t ZTF

python coverage.py -s ../data/GW190425/bayestar.fits.gz -f ../data/GW190425/Gattini_fields_bayestar.dat -l ../data/GW190425/transients_Gattini.dat -o ../output/GW190425/Gattini_bayestar.pdf -t Gattini

python coverage.py -s ../data/GW190425/LALInference.fits.gz -f ../data/GW190425/Gattini_fields.dat -l ../data/GW190425/transients_Gattini.dat -o ../output/GW190425/Gattini_LALInf.pdf -t Gattini

