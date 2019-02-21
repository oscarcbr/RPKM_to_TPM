# RPKM_to_TPM
Files to convert RPKMs to TPMs (see Wagner et al. 2012)

This repository contains a python (v2.7) script to convert RPKMs obtained from 
rpkmforgenes to an equivalent file with TPM values calculated as defined by
Wagner et al. 2012. 

How to execute?

python RPKM_to_TPM_rpkmforgenes.py <input_file_from_rpkmforgenes>

Input:

"input_file_from_rpkmforgenes" is a text file with rpkm values obtained from rpkmforgenes for various samples.

Output:

"input_file_from_rpkmforgenes"+.tpm is the output file in the same format with TPM values.

Example:

python RPKM_to_TPM_rpkmforgenes.py ensembl_rpkms_DEC13A.txt

