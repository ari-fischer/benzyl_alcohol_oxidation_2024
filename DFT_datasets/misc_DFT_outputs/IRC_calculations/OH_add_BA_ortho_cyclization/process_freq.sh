### extract outputs from freq calc output
grep 'Frequency:' OUT.out > Freq_out.txt
cut -c 14-100 Freq_out.txt > Frequencies.txt

grep 'Enthalpy' OUT.out > H.txt
cut -c 30-40 H.txt > H_out.txt

grep 'Entropy' OUT.out > S.txt
cut -c 30-40 S.txt > S_out.txt

grep 'Eigenvalues --' OUT.out > inertia_out.txt
cut -c 22-200 inertia_out.txt > inertia.txt

grep 'Red. Mass:' OUT.out > mass_out.txt
cut -c 14-100 mass_out.txt > masses.txt