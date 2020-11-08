#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
July 15 2020
smalkaram@wvstateu.edu

* This script is mirrored from,
* https://github.com/sridharacharya/bioStructureTools/blob/master/Structure.py
* Check the latest version of this script at this repo

Function: 
  This script
    1. takes a pdb file or pdb id as input
    2. optionally adds hydrogens
    3. then writes user selected residue(s) including connect records to output file

Dependencies :
    openbabel (optional: for converting formats or adding hydrogens)
    PySimpleGUI (optional: if gui is used)
    biopython (optional: if pdb-ids are used)
'''


class Atom:
    def __init__(self, line = None):
        import re
        if line is not None:
            self.line = line
        else:
            self.line = None

    def parse(self,line = None, opt=None):
        import re
        if line is not None:
            self.line = line

        if self.line is not None:
            z = re.findall('^([HA][ET][TO][AM][T\s][M\s])(.....).(....)(.)(...).(.)(....)(.)...(........)(........)(........)(......)(......)......(....)*(..)*(..)*', self.line)

            # The fields matched are
            #atomtyp #atomno #resna #chainna #resno #x #y #z #occupancy #tempfac #segmentid #elementSymbol #charge

            if len(z) == 0:
                return None
            elif opt is None:
                return z[0]
            elif opt == "coord":
                return [float(x) for x in z[0][8:11]]
            elif opt == "type":
                return z[0][0]
            elif opt == "atomno":
                return int(z[0][1])
            elif opt == "atomna":
                return z[0][2]
            elif opt == "altloc":
                return z[0][3]
            elif opt == "resna":
                return z[0][4]
            elif opt == "chainna":
                return z[0][5]
            elif opt == "resno":
                return int(z[0][6])
            elif opt == "inscode":
                return z[0][7]
            elif opt == "x":
                return float(z[0][8])
            elif opt == "y":
                return float(z[0][9])
            elif opt == "z":
                return float(z[0][10])
            elif opt == "occupancy":
                return z[0][11]
            elif opt == "tempfac":
                return z[0][12]
            elif opt == "segmentid":
                return z[0][13]
            elif opt == "elementSymbol":
                return z[0][14]
            elif opt == "charge":
                return float(z[0][15])
        else:
            return None



    def coord(self, line = None):
        if line is not None:
            self.line = line
        return self.parse(opt="coord")
    def type(self, line = None):
        if line is not None:
            self.line = line
        return self.parse(opt="type")
    def atomno(self, line = None):
        if line is not None:
            self.line = line
        return self.parse(opt="atomno")
    def resna(self, line = None):
        if line is not None:
            self.line = line
        return self.parse(opt="resna")
    def resno(self, line = None):
        if line is not None:
            self.line = line
        return self.parse(opt="resno")
    def chainna(self, line = None):
        if line is not None:
            self.line = line
        return self.parse(opt="chainna")
    def x(self, line = None):
        if line is not None:
            self.line = line
        return self.parse(opt="x")
    def y(self, line = None):
        if line is not None:
            self.line = line
        return self.parse(opt="y")
    def z(self, line = None):
        if line is not None:
            self.line = line
        return self.parse(opt="z")
    def occupancy(self, line = None):
        if line is not None:
            self.line = line
        return self.parse(opt="occupancy")
    def tempfac(self, line = None):
        if line is not None:
            self.line = line
        return self.parse(opt="tempfac")
    def segmentid(self, line = None):
        if line is not None:
            self.line = line
        return self.parse(opt="segmentid")
    def elementSymbol(self, line = None):
        if line is not None:
            self.line = line
        return self.parse(opt="elementSymbol")
    def charge(self, line = None):
        if line is not None:
            self.line = line
        return self.parse(opt="charge")
        
class Molecule(Atom):
    def __init__(self, FILE=None, ID=None, H=False):
       '''
       Pass optionally the pdb ID or PDB file
       and read the file
       '''
       self.id = ID
       self.file = FILE
       self.H = H
       self.selLines = []
       self.selAtoms = []
       self.selResidues = []
       self.selChains = []
       self.selConectLines= []
       self.selAtomLines= []
       self.ar = Atom()

       from datetime import datetime
       self.exectime = str(datetime.now())

       if ID is not None:
           FILE = self.getPDB(ID)
           self.file = FILE
       if FILE is not None:
           self.readPDB(FILE)

    def getPDB(self, ID=None):
        '''
        Retrives a PDB file from RCSB when ID is supplied
        or OBJECT.id is defined
        '''
        from Bio.PDB import PDBList
        if ID is not None :
            return PDBList().retrieve_pdb_file(ID, pdir = '.', file_format = 'pdb')
        elif self.id is not None:
            return PDBList().retrieve_pdb_file(self.id, pdir = '.', file_format = 'pdb')

    def addH(self, INPUT=None, OUTPUT=None, ID=None, oFORMAT="pdb", TITLE="None"):
        '''
        run openbabel to add hydrogens
        inputs: 
            INPUT=<string> (filename)
                (or) OBJECT.file defined
            OUTPUT=<string> (filename)
                (or) OBJECT.output defined
            ID=<string> (filename)
                (or) OBJECT.id defined
        outputs:
            PDB (OUTPUT) file with hydrogens added
        requirements:
            module openbabel
        '''
        if ID is not None:
            self.id = ID
        if INPUT is not None:
            self.file = INPUT
        if OUTPUT is not None:
            self.output = OUTPUT
        if self.id is None:
            print("Add .id first")
            return 0
        if self.file is not None:
            from openbabel import openbabel
            import re
            obMol = openbabel.OBMol()
            obConv = openbabel.OBConversion()

            if oFORMAT == "mol2":
                obConv.SetInAndOutFormats("pdb", "mol2")
            else :
                obConv.SetInAndOutFormats("pdb", "pdb")
            obConv.ReadFile(obMol, self.file)
            obMol.AddHydrogens()
            if TITLE is not None:
                obMol.SetTitle(TITLE)
            if self.output is not None:
                outMDL = obConv.WriteFile(obMol, self.output)
            else:
                outMDL = obConv.WriteFile(obMol, self.id + "_h.pdb")
                self.output = self.id + "_h.pdb"
            if outMDL is not None:
                self.H = True
        else:
            print("File not attached")
            print("use .file to add file name")
            print("or use .readPDB(<FILE>) to attach file")

    def readPDB(self, FILE=None):
        import re
        """
        Parse PDB file
        """
        #https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
        LINES={}
        CHAINS = {}
        RESNO={}
        RESNOEND={}
        RESNA={}
        RMAP={}  # NEW: to map the original residue and chain number to gromacs modified numbers
        ATOMNO={}
        HETATM=set()
        CONECT={}
        MASTER=0
        VALID = [0,0]
        natoms = 0
        nconnect = 0
        firstAtomno = 99999999
        lastAtomno = 0 

        if FILE is None:
            if self.file is not None:
                FILE = self.file
            else:
                print("Pass a FILE=<FILE> argument")
                return 0

        with open(FILE) as f:
            lineno=0

            ATOMREC = re.compile('^([AH][ET][TO][AM][T\s][M\s])(.....)......(...).(.)(....)')
            CONEREC = re.compile('^CONECT')
            COMPREC = re.compile('^COMPND\s+\d+MOLECULE:\s+(.*)')
            RESOREC = re.compile('^REMARK\s+\d+RESOLUTION\.\s+(.*)')
            MASTREC = re.compile('^MASTER')

            SERIAL_RESNO = 0
            LAST_RESNO = -1

            for line in f:
                line = line.replace('\n', '')   # Careful, only use this method, to just remove newline. Otherwise there will be parsing problems
                lineno += 1

                # Read ATOM or HETATM records
                #z1 = re.match('^([AH][ET][TO][AM][T\s][M\s])(.....)......(...).(.)(....)', line)
                z1 = ATOMREC.match(line)
                z2 = None
                # First look for CONECT record
                #if re.match('^CONECT', line):
                if CONEREC.match(line):
                    # Then check all matches for connections
                    z2 = re.findall(r'(.{5})', line.replace("CONECT", ""))
                #z3 = re.match('^MASTER', line)
                z3 = MASTREC.match(line)
                z4 = COMPREC.match(line)
                z5 = RESOREC.match(line)

		# ATOM or HETATM record
                if z1 is not None:
                    natoms += 1
                    LINES[lineno] = line
                    atomtyp = z1[1].strip()
                    atomno = int(z1[2].strip())
                    resna = z1[3].strip()
                    #chainna = z1[4].strip()
                    chainna = z1[4] # This takes care of empty chain names
                    #resno = int(z1[5].strip())
                    resno = chainna + z1[5].strip()

                    this_resno = int(z1[5].strip())

                    #if LAST_RESNO == 0 :
                    #    LAST_RESNO = this_resno

                    if SERIAL_RESNO == 0:
                        SERIAL_RESNO = this_resno
                        self.StartResidue = this_resno
                    elif this_resno != LAST_RESNO:
                        SERIAL_RESNO += 1
                        #RMAP[resno] = ' ' + str(this_resno - FACTOR + 1)

                    RMAP[resno] = ' ' + str(SERIAL_RESNO)
                    LAST_RESNO = this_resno 
                        
                    if atomno < firstAtomno : 
                        firstAtomno = atomno
                    if atomno > lastAtomno :
                        lastAtomno = atomno

                    # Store chain name and residue numbers
                    if chainna not in CHAINS:
                        CHAINS[chainna] = set() # Create if not present
                    CHAINS[chainna].add(resno)

                    '''
                    # Store an index for atomno 
                    ANDX = atomno // 100 + 1
                    if ANDX not in ATOMNO :
                        ATOMNO[ANDX] =  lineno
                    # Instead of index as above store actual atom number : lines
                    # as below
                    '''
                    ATOMNO[atomno] =  lineno

		    # If HETATM, add to the list
                    if atomtyp == "HETATM":
                        if resno not in HETATM:
                            HETATM.add(resno)

                    # This marks the beginning of residue
                    if resno not in RESNO :
                        RESNO[resno] = lineno
                        RESNA[resno] = resna

                    # This marks the end of residue
                    RESNOEND[resno] = lineno

                elif z2 is not None:
                    # Read CONECT records
                    nconnect += 1
                    for ITEM in z2:   # Number of atoms in CONECT record
                        if re.match(r'^\s*$', ITEM):
                            continue
                        ITEM = int(ITEM)
                        if not ITEM in CONECT:
                            CONECT[ITEM] = set()   # Create a set for this atom
                        CONECT[ITEM].add(lineno) # Add all lines # related to this atom
                    LINES[lineno] = line

                elif z3 is not None:
                    # Read MASTER
                    LINES[lineno] = line
                    MASTER = lineno

        # Sort residues in CHAINS
        # Note the sets are converted to lists here
        RC = re.compile(r'^(.)(\d+)')
        for C in CHAINS.keys():
            CHAIN = list(CHAINS[C])
            RNO = [ int(RC.findall(R)[0][1]) for R in CHAIN]
            CHAINS[C] = [CHAIN[RNO.index(R)] for R in sorted(RNO)]

        self.ATOMNO = ATOMNO
        self.RESNO = RESNO
        self.RESNOEND = RESNOEND
        self.RESNA = RESNA
        self.LINES = LINES 
        self.CHAINS = CHAINS
        self.PROTEIN = [I for I in self.RESNO.keys() if I not in HETATM]
        self.HETATM = HETATM
        self.CONECT = CONECT
        self.MASTER = MASTER
        self.nAtoms = natoms
        self.nConnect = nconnect
        self.firstAtomno = firstAtomno
        self.lastAtomno = lastAtomno
        self.RMAP = RMAP

        # Check if atleast some ATOM/HETATM or CONECT records are found
        if natoms < 1 :
          print("No ATOM lines found")
        if nconnect < 1 :
          print("No CONECT lines found")

    def summary(self):
        print("PDB ID = " + str(self.id))
        print("Number of chains = " + str(self.nChains()))
        print("    " + ",".join(sorted(self.CHAINS.keys())))
        print("Number of residues = " + str(self.nResidues()))
        print("    " + ",".join([ str(I) + ":" + str(len(self.CHAINS[I]))  for I in sorted(self.CHAINS.keys())]))
        print("Number of heteroatom residues = " + str(len(self.HETATM)))
        print("    " + ",".join([ str(I) + ":" + self.RESNA[I]  for I in sorted(self.HETATM) ]))
        print("Number of atoms = " + str(self.nAtoms))
        print("Number of connect records = " + str(self.nConnect))

    def lRange(self, L,R):
        D = {'A': 65, 'B': 66, 'C': 67, 'D': 68, 'E': 69, 'F': 70, 'G': 71, 'H': 72, 'I': 73, 'J': 74, 'K': 75, 'L': 76, 'M': 77, 'N': 78, 'O': 79, 'P': 80, 'Q': 81, 'R': 82, 'S': 83, 'T': 84, 'U': 85, 'V': 86, 'W': 87, 'X': 88, 'Y': 89, 'Z': 90}
        return list(map(chr, range(int(D[L]),int(D[R])+1)))

    def checkIfResiduePresent(self, R):
        if str(R) in self.RESNO:
            return 1
        else:
            print("Warning: Residue " + str(R) + " Not Found")
            return 0

    def checkIfChainPresent(self, C):
        if str(C) in self.CHAINS:
            return 1
        else:
            print("Warning: Chain " + str(C) + " Not Found")
            return 0


    def resna2no(self, resna = None):
        """
        Get the residue numbers for a specific name
        """
        if resna is None:
            print("Supply a residue name as argument")
        else:
            return [k for k,v in self.RESNA.items() if v == resna]

    def nChains(self):
        """
        Get the number of chains
        """
        if self.CHAINS is not None:
            return len(self.CHAINS.keys())

    def nResidues(self):
        """
        Get the number of residue
        """
        if self.RESNO is not None:
            return len(self.RESNO.keys())

    def Chains(self):
        """
        Print chains
        """
        if self.CHAINS is not None:
            return self.CHAINS.keys()

    def Resnos(self):
        """
        Print residue numbers
        """
        if self.RESNO is not None:
            return self.RESNO.keys()

    def Hetatms(self):
        """
        Print heteroatom residues
        """
        if self.HETATM is not None:
            return self.HETATM

    def ResInChain(self, CHAIN=None):
        """
        Print residue numbers in a chain
        """
        if self.CHAINS is not None:
            return self.CHAINS[CHAIN]

    def whichChain(self, RESIDUE=None):
        """
        Print the chain of a particular residue
        """
        if self.RESNO is not None:
            return self.ar(line = self.RESNO[RESIDUE], opt = "chainna")

    def whichResidue(self, RESIDUE=None):
        """
        Print residue name for a particular number
        """
        if self.RESNO is not None:
            return self.RESNA[RESIDUE]
    
    def getConnectRecords(self, ATMNOS=None):
        """
        get connect records
        """
        self.selConectLines = []
        if ATMNOS is not None:
            self.selConectLines = sorted({J for I in ATMNOS if I in self.CONECT for J in self.CONECT[I]})
        elif self.selAtoms is not None:
            self.selConectLines = sorted({J for I in self.selAtoms if I in self.CONECT for J in self.CONECT[I]})
        #return self.selConectLines

    def selectDialog(self):
        import PySimpleGUI as sg

        # This is overall layout
        LAYOUT = []

        CHAINS = self.CHAINS.keys()
        CHAINS_VISIBLE = dict(zip(CHAINS, [0] * len(CHAINS)))
        HETATM_VISIBLE = 0

        IS_FIRST = 0

        BUTTONS = []

        # For each Chain
        for C in CHAINS :
            RES = [str(I) + ":" + str(self.RESNA[I]) for I in self.CHAINS[C]]
            layout = []
            BUTTONS.append(sg.Button(C))
            layout.append([sg.Text("chain " + C)])
            layout.append([sg.Listbox(values=["None", "All", "Selected"], select_mode='single', size=(10,3))])
            layout.append([sg.Text('Residues')])
            layout.append([sg.Listbox(values= RES, select_mode='multiple', size=(10,10))])
            layout.append([sg.Checkbox('Invert', default=False)])

            if IS_FIRST == 0:
                LAYOUT.append(sg.Column(layout, visible=True, key="chain-" + C))
                CHAINS_VISIBLE[C] = 1 
                IS_FIRST = 1
            else:
                LAYOUT.append(sg.Column(layout, visible=False, key="chain-" + C))

        # For heteroatoms
        HETERO  = [str(I) + ":" + str(self.RESNA[I]) for I in self.HETATM]
        layout = []
        BUTTONS.append(sg.Button("HET"))
        layout.append([sg.Text("Heteroatoms")])
        layout.append([sg.Listbox(values=["None", "All", "Selected"], select_mode='single', size=(10,3))])
        layout.append([sg.Text('Residues')])
        layout.append([sg.Listbox(values=HETERO, select_mode='multiple', size=(10,10))])
        layout.append([sg.Checkbox('Invert', default=False)])

        LAYOUT.append(sg.Column(layout, visible = False, key="HETATM"))

        # For submit
        BUTTONS2 = [sg.Submit(), sg.Cancel()]

        layout = [sg.Text("Click on these categories to hide/unhide them")]


        LAYOUT2 = [layout, BUTTONS, LAYOUT, BUTTONS2] 

        # window
        window = sg.Window("Select Residues").Layout(LAYOUT2)

        # check response
        while True:
            button, values = window.Read()
            if button in CHAINS :
                if CHAINS_VISIBLE[button] == 0 :
                    CHAINS_VISIBLE[button] = 1 
                    window[f'chain-{button}'].update(visible=True)
                else :
                    CHAINS_VISIBLE[button] = 0
                    window[f'chain-{button}'].update(visible=False)
            elif button == "HET" :
                if HETATM_VISIBLE == 0:
                    HETATM_VISIBLE = 1
                    window[f'HETATM'].update(visible=True)
                else :
                    HETATM_VISIBLE = 0
                    window[f'HETATM'].update(visible=False)
            else :
                if button != "Submit" :
                    values = dict(zip(range(0, (self.nChains() + 1)* 3), [[],[],False] * (self.nChains() + 1)))
                    window.Close()
                    break
                if button is not None :
                    window.Close()
                    break
                if button == WIN_CLOSED :
                    window.Close()
                    break
        return values

    def selectByDialog(self, TEXT="Select residues"):

        if self.PROTEIN is None and self.HETATM is None:
            print("Empty object")
            return None

        SELECTION = list(self.selectDialog().values())
        CHAINS = list(self.CHAINS.keys())
        rPROTEIN = []

        # Return from selectDialog has 3 elements/chain and heteroatoms, so run a loop to capture these serially
        # For chains
        for C in range(0, len(CHAINS)):

            START = C * 3
            END = START + 3

            TYPE, SEL, INVERT = SELECTION[START:END]

            if len(TYPE) == 0:
                TYPE = None
                INVERT = False
                SEL = []
            else :
                TYPE = TYPE[0]

            CHAIN = CHAINS[C]
            PROTEIN = []

            if INVERT :
                if TYPE == "None":
                    #PROTEIN = list(self.CHAINS[CHAIN])
                    PROTEIN = self.CHAINS[CHAIN]
                elif TYPE == "All":
                    PROTEIN = []
                elif TYPE == "Selected":
                    TEMP = [K.split(':')[0] for K in SEL]
                    PROTEIN = [ L for L in self.PROTEIN if L not in TEMP ]
            else:
                if TYPE == "None":
                    PROTEIN = []
                elif TYPE == "All":
                    #PROTEIN = list(self.CHAINS[CHAIN])
                    PROTEIN = self.CHAINS[CHAIN]
                elif TYPE == "Selected":
                    PROTEIN = [K.split(':')[0] for K in SEL]

            rPROTEIN = rPROTEIN + PROTEIN

        # For heteroatoms
        START = C * 3
        END = START + 3
        TYPE, SEL, INVERT = SELECTION[START:END]
        HETERO = []

        if INVERT :
            if TYPE == "None":
                HETERO = list(self.HETATM)
            elif TYPE == "All":
                HETERO = []
            elif TYPE == "Selected":
                TEMP = [K.split(':')[0] for K in SEL ]
                HETERO = [ L for L in self.HETATM if L not in TEMP ]
        else:
            if TYPE == "None":
                HETERO = []
            elif TYPE == "All":
                HETERO = list(self.HETATM)
            elif TYPE == "Selected":
                HETERO = [K.split(':')[0] for K in SEL]

        self.selResidues = rPROTEIN + HETERO
        return rPROTEIN + HETERO

    def selectByText(self, TEXT = None):
        '''
        Provide a TEXT string to make a selection with following methods
        'A'  = selects chain A
        'A-G' = selects chain A to G
        'A123' = selects residue 123 in chain A
        '10-20' = selects residues 10 to 20 in all chains
        'A123,B126' = selects residue 123 in chain A and residue 126 in chain B
        Any other combinations, where selection elements are separated by comma
        '''
        import re
        SELECTION = [] 
        CHAINS = self.CHAINS.keys()
        RANGE = re.compile('[^\:]*:[^\:]*')

        CHAIN = re.compile(r'^[A-Z\s]$') # Chain identifier can be empty
        NUMBER = re.compile(r'^\d+$')
        RESIDUE = re.compile(r'^[A-Z\s]\d+$', re.IGNORECASE) # Chain identifier can be empty


        # split each term
        for SEL in TEXT.split(","):
            # If this is a range
            if RANGE.match(SEL) :
                L,R = SEL.split(":")
                # if this is a number
                if NUMBER.match(L):
                    for C in CHAINS :
                        for R in range(int(L), int(R) +1):
                            self.checkIfResiduePresent(C + str(R)) and SELECTION.append(C + str(R))
                # if this is a chain
                elif CHAIN.match(L):
                    for C in self.lRange(L,R):
                        if self.checkIfChainPresent(C) :
                            SELECTION = SELECTION + list(self.CHAINS[C])
                # if this is a residue
                elif RESIDUE.match(L):
                    Lc = re.findall("^[A-Z]", L)[0] 
                    if Lc is not re.findall("^[A-Z]", R)[0]:
                        print("Residue selection: The two chains in the range must be same")
                        exit(1)
                        
                    for R in range(int(re.sub(r"^[A-Z]","", L)), int(re.sub(r"^[A-Z]", "", R)) + 1):
                        self.checkIfResiduePresent(Lc + str(R)) and SELECTION.append(Lc + str(R))

            # If this is a single
            else:
                # if this is a number
                if NUMBER.match(SEL):
                    for C in CHAINS :
                        self.checkIfResiduePresent(C + str(SEL)) and SELECTION.append(C + str(SEL))
                # if this is a chain
                elif CHAIN.match(SEL):
                    if self.checkIfChainPresent(SEL) :
                        SELECTION = SELECTION + list(self.CHAINS[SEL])
                # if this is a residue
                elif RESIDUE.match(SEL):
                    self.checkIfResiduePresent(SEL) and SELECTION.append(SEL)
        self.selResidues = SELECTION
        return SELECTION

    def selectionExclude(self, INCLUDE = None, EXCLUDE = None):
        assert type(INCLUDE) is list
        assert type(EXCLUDE) is list
        return list(set(INCLUDE) - set(EXCLUDE))

    def changeResno(self, RE_RESNO, line = None, newResno = None):
        '''
        This function takes atom record and changes the residue numbers 
        according to a residue number supplied
        '''
        import re
        assert line is not None, "line should be supplied"

        #ATOM     90  ND2 ASN A  95     -28.481  73.958  89.395  1.00 16.41           N
        z1 = RE_RESNO.findall(line)
        RE_RESNO = re.compile('^([AH][ET][TO][AM][T\s][M\s]...............)(.....)(.*)')
        return "{0:<21s}{1:>5s}{2:<s}". format(z1[0][0], newResno, z1[0][2])


    def changeResna(self, RE_RESNAME, line = None, newResname = None):
        '''
        This function takes atom record and changes the residue name
        according to a residue name supplied
        '''
        import re
        assert line is not None, "line should be supplied"

        #ATOM     90  ND2 ASN A  95     -28.481  73.958  89.395  1.00 16.41           N
        z1 = RE_RESNAME.findall(line)
        return "{0:<17s}{1:<3s}{2:<s}". format(z1[0][0], newResname, z1[0][2])

    def changeChainna(self, RE_CHAINNAME, line = None, refMol= None, newChainna=None):
        '''
        This function takes atom record and changes the atom numbers 
        according to chain name supplied
        '''
        import re
        assert line is not None, "line should be supplied"

        assert refMol is not None or newChainna is not None, "refMol or newChainna should be supplied"

        #ATOM     90  ND2 ASN A  95     -28.481  73.958  89.395  1.00 16.41           N
        z1 = RE_CHAINNAME.findall(line)

        if newChainna is not None :
        #newChainname = refMol.chainna(refMol.LINES[refMol.ATOMNO[self.atomno(line)]])
            return "{0:<21s}{1:<1s}{2:<s}". format(z1[0][0], 
            newChainna,
            z1[0][2])
        else :
            return "{0:<21s}{1:<1s}{2:<s}". format(z1[0][0], 
            refMol.chainna(refMol.LINES[refMol.ATOMNO[self.atomno(line)]]),
            z1[0][2])


    def changeAtomnoInAtomRecord(self, line = None, lastAtomno = None):
        '''
        This function takes atom record and changes the atom numbers 
        according to a starting atom number supplied
        inputs:
            line = connect record
            lastAtomno = last atom number from which new atom numbers should start
        outputs:
            newline = changed atom record
        '''
        import re
        if line is not None:
            self.line = line
        z = re.findall('^([HA][ET][TO][AM][T\s][M\s])(.....)(.*)', self.line)

        if lastAtomno is not None:
            self.lastAtomno = lastAtomno

        newatomno = int(z[0][1]) + self.lastAtomno
        return "{0:<6s}{1:>5d}{2:<s}". format(z[0][0], newatomno, z[0][2])

    def changeAtomnoInConectRecord(self, line = None, lastAtomno = None):
        '''
        This function takes a connect record and changes the atom numbers 
        according to a starting atom number supplied
        inputs:
            line = connect record
            lastAtomno = last atom number from which new atom numbers should start
        outputs:
            newline = changed connect record
        '''
        import re
        if line is not None:
            self.line = line

        if not re.match('^CONECT', self.line):
            print("This is not a connect record")
            exit(2)

        if lastAtomno is not None:
            self.lastAtomno = lastAtomno

        # Then check all matches for connections
        z2 = re.findall(r'(.{5})', self.line.replace("CONECT", ""))

        newline = "CONECT"

        for atomno in z2 :
            #print("atomno = " + str(atomno))
            newline = newline + str("{0:>5d}". format(int(atomno) + self.lastAtomno))

        return newline


    def selectResidues(self, resno = None):
        """
        Prints the selected residue number including the connect records
        """
        import re
        if resno is not None:
            self.selResidues = resno
        assert self.selResidues is not None, "Needs either residue numbers list as argument or .selResidues defined"

        # Get a list of line numbers for these residues 
        self.selAtomLines = sorted({
            I
            for RESNO in self.selResidues
            for I in range(self.RESNO[RESNO], self.RESNOEND[RESNO]+1) 
            # Check if these lines are indeed for this residue selection
            if str(self.ar.resno(self.LINES[I])) == re.sub("^.", "", RESNO)
            })
            #if self.ar.resno(self.LINES[I]) == RESNO

        # Now select the atomno's for the above residue list
        self.selAtoms = sorted({
            self.ar.atomno(self.LINES[I]) 
            for I in self.selAtomLines
            })

        # Get connect records
        self.getConnectRecords()

        # Create selected lines based on atom records and connect records
        self.selLines =  sorted(set(self.selAtomLines + self.selConectLines))

        #return self.selLines

    def write(self, FILE="output.pdb", changeAtomnos = False, lastAtomno = None, changeResna=False, newResname=None, changeResno=False, changeChainna=False, newResno=None, refPDB=None, newChainna=None):
        """
        write selected residues to file
        Note: Even though there are several ways to renumber residues, and rename residue names and chain names, the way gromacs ignores these original identities makes it very complex to properly usethese methods. So it is advisable to just use the residue mappings rather than attempting to renumber or rename.
        """
        if self.id is None:
            self.id = "UNK"
        if self.selLines is not None:
            import re
            with open(FILE, "w") as fp:
                fp.write('HEADER     Created using PDB.Structure script on : ' + self.exectime + '\n')
                fp.write('TITLE     From PDB ID: ' + self.id + '\n')
                '''
                for I in re.findall('(.{,66})', ' '.join([str(I) + ":" + self.RESNA[I] for I in self.selResidues])):
                    if len(I) > 1:
                        fp.write('COMPND ' + "     " + I + '\n')
                '''
                if changeResna:
                    assert changeAtomnos is False, "Cant do both change Atomnos and Resna together"
                    assert changeResno is False, "Cant do both change Resno and Resna together"
                    assert changeChainna is False, "Cant do both change Resna and Channa together"
                    assert newResname is not None, "Must supply newResname"
                    RE_RESNAME = re.compile('^([AH][ET][TO][AM][T\s][M\s]...........)(...)(.......*)')
                if changeChainna:
                    assert changeAtomnos is False, "Cant do both change Atomnos and Channa together"
                    assert changeResno is False, "Cant do both change Resno and Channa together"
                    assert changeResna is False, "Cant do both change Channa and Resna together"
                    assert refPDB is not None or newChainna is not None, "Must supply a refPDB or newChainna"
                    if refPDB is not None:
                        refMol = Structure.Molecule(FILE=refPDB)
                    RE_CHAINNAME = re.compile('^([AH][ET][TO][AM][T\s][M\s]...............)(.)(.*)')

                if changeResno:
                    assert changeAtomnos is False, "Cant do both change Atomnos and Resno together"
                    assert changeResna is False, "Cant do both change Resno and Resna together"
                    assert changeChainna is False, "Cant do both change Resno and Chainna together"
                    assert changeResno is not None, "Must supply changeResno"
                    RE_RESNO = re.compile('^([AH][ET][TO][AM][T\s][M\s]...............)(.....)(.*)')

                if changeAtomnos :
                    if lastAtomno is None or int(lastAtomno) <= 0 :
                        if self.lastAtomno is None :
                            print("argument lastAtomno should be supplied and should be  >= 1")
                            return 0
                        else :
                            lastAtomno = self.lastAtomno
                    for I in self.selAtomLines:
                        fp.write( self.changeAtomnoInAtomRecord(line = self.LINES[I], lastAtomno = lastAtomno) + '\n')

                    for I in self.selConectLines:
                        fp.write( self.changeAtomnoInConectRecord(line = self.LINES[I], lastAtomno = lastAtomno) + '\n')
                else :
                    #for I in self.selAtomLines:
                    #    fp.write(self.LINES[I] + '\n')

                    #for I in self.selConectLines:
                    #    fp.write(self.LINES[I] + '\n')

                    # Write selected lines
                    if changeResna :
                        for I in self.selLines:
                            fp.write(self.changeResna(RE_RESNAME, self.LINES[I], newResname) + '\n')
                    elif changeResno :
                        for I in self.selLines:
                            fp.write(self.changeResno(RE_RESNO, self.LINES[I], newResno) + '\n')
                    elif changeChainna :
                        if newChainna is not None :
                            for I in self.selLines:
                                fp.write(self.changeChainna(RE_CHAINNAME, line=self.LINES[I], newChainna=newChainna) + '\n')
                        elif refMol is not None :
                            for I in self.selLines:
                                fp.write(self.changeChainna(RE_CHAINNAME, line=self.LINES[I], refMol=refMol) + '\n')
                    else :
                        for I in self.selLines:
                            fp.write(self.LINES[I] + '\n')
                fp.write('END\n')
        else:
            print("Object has no data to write")


class Complex(Molecule):
    def __init__(self, proteinfile=None, ligandfile=None):
       '''
       Pass optionally the pdb ID or PDB file
       and read the file
       '''
       self.proteinfile = proteinfile
       self.ligandfile = ligandfile
       self.protein = None
       self.ligand = None
       self.read(proteinfile = self.proteinfile, ligandfile = self.ligandfile)

    def read(self, proteinfile = None, ligandfile = None):
        if proteinfile is not None:
            self.proteinfile = proteinfile
        if ligandfile is not None:
            self.ligandfile = ligandfile
        if self.proteinfile is not None :
            self.protein = Molecule(FILE=self.proteinfile)
        if self.ligandfile is not None:
            self.ligand = Molecule(FILE=self.ligandfile)

    def combine(self, proteinsel = None, ligandsel = None):
        if proteinsel is not None:
            self.protein.selResidues = proteinsel
        if ligandsel is not None:
            self.ligand.selResidues = ligandsel
        if self.protein.selResidues is None or self.protein.selResidues == []:
            self.protein.selResidues = list(self.protein.RESNO.keys())
        if self.ligand.selResidues is None or self.ligand.selResidues == []:
            self.ligand.selResidues = list(self.ligand.RESNO.keys())

        self.protein.selectResidues(resno = self.protein.selResidues)
        self.ligand.selectResidues(resno = self.ligand.selResidues)

    def write(self, FILE = "complex.pdb") :
        import re
        if FILE is not None:
            self.output = FILE

        with open(self.output, "w") as fp :
            for I in re.findall('(.{,66})', self.proteinfile + " + " + self.ligandfile):
                if len(I) > 1:
                    fp.write('TITLE  ' + "     " + I + '\n')
            '''
            for I in re.findall('(.{,66})', ' '.join([str(I) + ":" + self.protein.RESNA[I] for I in self.protein.selResidues])):
                if len(I) > 1:
                    fp.write('COMPND  ' + "     " + I + '\n')
            for I in re.findall('(.{,66})', ' '.join([str(I) + ":" + self.ligand.RESNA[I] for I in self.ligand.selResidues])):
                if len(I) > 1:
                    fp.write('COMPND  ' + "     " + I + '\n')
            '''
            for I in sorted(set(self.protein.selAtomLines)):
                fp.write(self.protein.LINES[I] + '\n')

            fp.write('TER\n')

            for I in sorted(set(self.ligand.selAtomLines)):
                fp.write(self.ligand.LINES[I] + '\n')

            for I in sorted(set(self.protein.selConectLines)):
                fp.write(self.protein.LINES[I] + '\n')
            for I in sorted(set(self.ligand.selConectLines)):
                fp.write(self.ligand.LINES[I] + '\n')

            fp.write('END\n')


def main(argv):
    try:
        opts, args = getopt.getopt(argv,":hs:if:r:x:o:")
    except getopt.GetoptError:
        print('Usage: Structure.py -i <pdb-id> or -f <pdb-file> -r <included resno> -x <excluded resno> -o <output-file-name>')
        print('    Use option -s for summary')
        sys.exit(2)
    FILE = None
    ID = None
    EXCLUDE = None
    SUMMARY = None
    RESIDUE = None
    OUTPUT = "output.pdb"
    for opt, arg in opts:
        if opt == '-h':
            print('Usage: Structure.py -i <pdb-id> or -f <pdb-file> -r <included resno> -x <excluded resno> -o <output-file-name>')
            print('    Use option -s for structure summary')
            print('    -h for this help')
            sys.exit()
        elif opt in ("-s"):
            SUMMARY = arg
        elif opt in ("-f"):
            FILE = arg
        elif opt in ("-i"):
            ID = arg
        elif opt in ("-x"):
            EXCLUDE = arg
        elif opt in ("-r"):
            RESIDUE = arg
        elif opt in ("-o"):
            OUTPUT = arg

    dat = Molecule(FILE=FILE, ID = ID)

    dat.readPDB(FILE=FILE)

    if SUMMARY is not None:
        dat.summary()
        exit(1)

    if RESIDUE is not None:
        SELECTED = dat.selectByText(RESIDUE)
    else:
        SELECTED = dat.selectByDialog()

    if EXCLUDE is not None:
        EXCLUDE = dat.selectByText(EXCLUDE)
        dat.selectResidues(resno = dat.selectionExclude(INCLUDE=SELECTED, EXCLUDE=EXCLUDE))
    else :
        dat.selectResidues(resno = SELECTED)

    dat.write(FILE=OUTPUT)

if __name__ == "__main__":
    import sys,getopt
    main(sys.argv[1:])
