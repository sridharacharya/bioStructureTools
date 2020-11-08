# bioStructureTools
This repository contains some utilities to parse/edit/write biomolecular structure files (mainly PDB files)

### Structure.py 
- Takes a pdb file or pdb id as input
- optionally adds hydrogens
- writes user selected residue(s) including connect records to output file

- Usage: Structure.py -i <pdb-id> or -f <pdb-file> -r <included resno> -x <excluded resno> -o <output-file-name>
  - Use option \-s for summary
