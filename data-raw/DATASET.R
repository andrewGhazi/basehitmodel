# The data in data-raw/ are compressed, barcode-anonymized count files. The files were
# compressed with:

# $ xz --version
# xz (XZ Utils) 5.2.5
# liblzma 5.2.5

# The data in demo/ are the first 1000 columns of the first three files i.e. uncompress the files with xz -kd and then:
# cut -f -1000 -d, Mapped_bcs_run1.csv > demo1.csv
