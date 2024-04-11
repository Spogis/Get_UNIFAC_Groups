UNIFAC_SMARTS = [
        ("CH3", "[CH3;X4]"),
        ("CH2", "[CH2;X4]"),
        ("CH", "[CH1;X4]"),
        ("C", ["[CH0;X4]", "[CH0;X3]"]),
        ("CH2=CH", "[CH2]=[CH]"),
        ("CH=CH", "[CH]=[CH]"),
        ("CH2=C", ["[CH2]=[C]", "[CH2]=[c]"]),
        ("CH=C", ["[CH]=[CH0]", "[CH]=[cH0]"]),
        ("ACH", "[cH]"),
        ("AC", "[cH0]"),
        ("ACCH3", "[c][CH3;X4]"),
        ("ACCH2", "[c][CH2;X4]"),
        ("ACCH", "[c][CH;X4]"),
        ('OH', "[OH]"),
        ('CH3OH', "[CH3][OH]"),
        ('H2O', "[OH2]"),
        ('ACOH', "[c][OH]"),
        ("CH3CO", "[CH3][CH0]=O"),
        ("CH2CO", "[CH2][CH0]=O"),
        ("CH=O", "[CH]=O"),
        ("CH3COO", "[CH3]C(=O)[OH0]"),
        ("CH2COO", "[CH2]C(=O)[OH0]"),
        ("HCOO", "[CH](=O)[OH0]"),
        ("CH3O", "[CH3][OH0]"),
        ("CH2O", "[CH2][OH0]"),
        ("CHO", "[CH][OH0]"),
        ("THF", "[CH2;R][OH0]"),
        ("CH3NH2", "[CH3][NH2]"),
        ("CH2NH2", "[CH2][NH2]"),
        ("CHNH2", "[CH][NH2]"),
        ("CH3NH", "[CH3][NH]"),
        ("CH2NH", "[CH2][NH]"),
        ("CHNH", "[CH][NH]"),
        ("CH3N", ["[CH3][N]", "[CH3][n]"]),
        ("CH2N", "[CH2][N]"),
        ("ACNH2", "[c][NH2]"),
        ("C5H5N", "n1[cH][cH][cH][cH][cH]1"),
        ("C5H4N", ["n1[c][cH][cH][cH][cH]1",
                   "n1[cH][c][cH][cH][cH]1",
                   "n1[cH][cH][c][cH][cH]1"]),
        ('C5H3N', ["n1[c][c][cH][cH][cH]1",
                   "n1[c][cH][c][cH][cH]1",
                   "n1[c][cH][cH][c][cH]1",
                   "n1[c][cH][cH][cH][c]1",
                   "n1[cH][c][c][cH][cH]1",
                   "n1[cH][c][cH][c][cH]1"]),
        ("CH3CN", "[CH3]C#N"),
        ("CH2CN", "[CH2]C#N"),
        ("COOH", "C(=O)[OH]"),
        ("HCOOH", "[CH](=O)[OH]"),
        ("CH2Cl", "[CH2]Cl"),
        ("CHCl", "[CH]Cl"),
        ("CCl", "[CH0]Cl"),
        ("CH2Cl2", "[CH2](Cl)Cl"),
        ("CHCl2", "[CH](Cl)Cl"),
        ("CCl2", "C(Cl)Cl"),
        ("CHCl3", "[CH](Cl)(Cl)Cl"),
        ("CCl3", "C(Cl)(Cl)(Cl)"),
        ("CCl4", "C(Cl)(Cl)(Cl)(Cl)"),
        ("ACCl", "[c]Cl"),
        ("CH3NO2", "[CH3][N+](=O)[O-]"),
        ("CH2NO2", "[CH2][N+](=O)[O-]"),
        ("CHNO2", "[CH][N+](=O)[O-]"),
        ("ACNO2", "[c][N+](=O)[O-]"),
        ("CS2", "C(=S)=S"),
        ("CH3SH", "[CH3][SH]"),
        ("CH2SH", "[CH2][SH]"),
        ("Furfural", "O=[CH]c1[cH][cH][cH]o1"),
        ("DOH", "[OH][CH2][CH2][OH]"),
        ("I", "[IH0]"),
        ("Br", "[BrH0]"),
        ("CH#C", "[CH]#C"),
        ("C#C", "C#C"),
        ("DMSO", "[CH3]S(=O)[CH3]"),
        ("ACRY", "[CH2]=[CH1][C]#N"),
        ("Cl(C=C)", "[$(Cl[C]=[C])]"),
        ("C=C", "[CH0]=[CH0]"),
        ("ACF", "[c]F"),
        ("DMF", "[CH](=O)N([CH3])[CH3]"),
        ("HCON(CH2)2", ["[CH](=O)N([CH2])[CH2]", "[CH](=O)N([CH2])[CH3]"]),
        ("CF3", "C(F)(F)F"),
        ("CF2", "C(F)F"),
        ("CF", "[C]F"),
        ("COO", ["[CH0](=O)[OH0]", "[cH0](=O)[oH0]"]),
        ("SiH3", "[SiH3]"),
        ("SiH2", "[SiH2]"),
        ("SiH", "[SiH]"),
        ("Si", "[Si]"),
        ("SiH2O", "[SiH2][OH0]"),
        ("SiHO", "[SiH][OH0]"),
        ("SiO", "[Si][OH0]"),
        ("NMP", "[CH3]N1[CH2][CH2][CH2]C(=O)1"),
        ("CCl3F", "C(Cl)(Cl)(Cl)F"),
        ("CCl2F", "C(Cl)(Cl)F"),
        ("HCCl2F", "[CH](Cl)(Cl)F"),
        ("HCClF", "[CH](Cl)F"),
        ("CClF2", "C(Cl)(F)F"),
        ("HCClF2", "[CH](Cl)(F)F"),
        ("CClF3", "C(Cl)(F)(F)F"),
        ("CCl2F2", "C(Cl)(Cl)(F)F"),
        ("CONH2", "C(=O)[NH2]"),
        ("CONHCH3", "C(=O)[NH][CH3]"),
        ("CONHCH2", "C(=O)[NH][CH2]"),
        ("CON(CH3)2", "C(=O)N([CH3])[CH3]"),
        ("CONCH3CH2", "C(=O)N([CH3])[CH2]"),
        ("CON(CH2)2", "C(=O)N([CH2])[CH2]"),
        ("C2H5O2", "[OH0;!$(OC=O);!R][CH2;!R][CH2;!R][OH]"),
        ("C2H4O2", ["[OH0;!$(OC=O);!R][CH;!R][CH2;!R][OH]",
                    "[OH0;!$(OC=O);!R][CH2;!R][CH;!R][OH]"]),
        ("CH3S", "[CH3]S"),
        ("CH2S", "[CH2]S"),
        ("CHS", "[CH]S"),
        ("MORPH", "[CH2]1[CH2][NH][CH2][CH2]O1"),
        ("C4H4S", "[cH]1[cH][s;X2][cH][cH]1"),
        ('C4H3S', ["[c]1[cH][s;X2][cH][cH]1",
                   "[cH]1[c][s;X2][cH][cH]1"]),
        ('C4H2S', ["[c]1[c][s;X2][cH][cH]1",
                   "[c]1[cH][s;X2][cH][c]1",
                   "[cH]1[c][s;X2][c][cH]1",
                   "[cH]1[c][s;X2][cH][c]1"]),
        ("NCO", "N=C=O"),
        ("H2COCH", ""),  # "[CH2]1[CH]O1"
        ("HCOCH", ""),  # "[CH]1[CH]O1"
        ("COCH", ""),  # "C1[CH]O1"
        ("H2COCH2", ""),  # "[CH2]1[CH2]O1"
        ("OCOCO", ""),  # "C(=O)OC(=O)"),
        ("(CH3O)2CO", ""),  # "[CH3]OC(=O)O[CH3]"
        ("(CH2O)2CO", ""),  # "[CH2]OC(=O)O[CH2]"
        ("CH3OCH2OCO", ""),  # "[CH3]OC(=O)O[CH2]"
        ("(CH2)2SU", "[CH2]S(=O)(=O)[CH2]"),
        ("CH2CHSU", "[CH2]S(=O)(=O)[CH]")]