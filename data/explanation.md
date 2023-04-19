## Data Folder

This folder contains data files for the radial profiles of objects through which our dark matter particle will propagate.

Each file contains a header row that firsts gives the mass of the object in grams, then the radius of the object in centimeters, then the masses in atomic mass units of each nuclear species whose composition is listed in the data file. Below the header, the first column lists the total mass fraction contained within each radial shell of the object. The second column lists the fraction of the total radius contained within each radial shell. The third column gives the temperature at each radius in Kelvin, and the fourth gives the density at each radius in grams per cubic centimeter. Each of the following columns gives the mass fraction composed by each of the listed nuclear species, at that radius. These species composition columns are listed in the same order as the masses of the species listed in the header.

The solar files "struct_b16_agss09_modified.dat" and "struct_b16_gs98_modified.dat" were modified from the files "struct_b16_agss09.dat" and "struct_b16_gs98.dat" retrieved from https://www.ice.csic.es/personal/aldos/Solar_Data.html. These files were modified to fit the format described in the above paragraph. To be precise, the large descriptive headers were removed and replaced with the one-line headers containing information about the total mass and radius of the Sun, as well as the masses of the listed species, as described above. Additionally, the columns for the radial pressure and luminosity were removed, as they are not necessary for this program.