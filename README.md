# bioStructureTools
This repository contains some utilities to parse/edit/write biomolecular structure files (mainly PDB files)

### Structure.py 
- Takes a pdb file or pdb id as input
- optionally adds hydrogens
- writes user selected residue(s) including connect records to output file
- The script is working but not completely error free.
- Use at your own risk, with adequate checking of output

```
- Usage: Structure.py -i <pdb-id> or -f <pdb-file> -r <included resno> -x <excluded resno> -o <output-file-name>
  - Use option \-s for summary
```
#### Dependencies :
- openbabel (optional: for converting formats or adding hydrogens)
- PySimpleGUI (optional: if gui is used)
- biopython (optional: if pdb-ids are used)
