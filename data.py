# -*- coding: utf-8 -*-

from rdkit.Chem import AllChem


def to_fingerprint(smiles):
    molecule = AllChem.MolFromSmiles(smiles)
    bit_vector = AllChem.GetMorganFingerprintAsBitVect(molecule, 4, nBits=512)
    fingerprint = bit_vector.ToBitString()

    return fingerprint


n = 15
m = 22
side_effects = ['нейтропения', 'гранулоцитопения', 'тромбоцитопения', 'гипопротромбинемия', 'тромбоцитопеническая пурпура', 'эозинофилия', 'транзиторная анемия', 'тахикардия', 'анемия', 'панцитопения', 'нарушения коагуляции', 'тромбоцитоз', 'тромбофлебит', 'агранулоцитоз', 'лейкопения', 'томбоцитопения', 'гемолитическая анемия', 'тромбоцитемия', 'удлинение ПВ', 'угнетение костного мозга', 'аритмия', 'супрессия костного мозга']
compounds = [['CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O', 'тахикардия', 'транзиторная анемия', 'тромбоцитопеническая пурпура', 'эозинофилия', 'лейкопения', 'нейтропения', 'агранулоцитоз'],
             ['CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O', 'лейкопения', 'нейтропения', 'тромбоцитопения', 'агранулоцитоз', 'анемия'],
             ['CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O.CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O.O.O.O.O.c1ccc(CNCCNCc2ccccc2)cc1', 'анемия', 'тромбоцитопения', 'лейкопения', 'нарушения коагуляции'],
             ['CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O'],
             ['Cc1onc(-c2ccccc2)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12', 'агранулоцитоз', 'нейтропения', 'угнетение костного мозга'],
             ['Cc1nnc(SCC2=C(C(=O)O)N3C(=O)[C@@H](NC(=O)Cn4cnnn4)[C@H]3SC2)s1', 'нейтропения', 'лейкопения', 'тромбоцитопения', 'тромбоцитемия'],
             ['CC1=C(C(=O)O)N2C(=O)[C@@H](NC(=O)[C@H](N)c3ccccc3)[C@H]2SC1', 'лейкопения', 'нейтропения', 'тромбоцитопения', 'удлинение ПВ', 'эозинофилия'],
             ['CC(=O)OCC1=C(C(=O)O)N2C(=O)[C@@H](NC(=O)Cc3cccs3)[C@H]2SC1'],
             ['Cn1nnnc1SCC1=C(C(=O)O)N2C(=O)[C@@H](NC(=O)[C@H](O)c3ccccc3)[C@H]2SC1', 'лейкопения', 'тромбоцитопения', 'нейтропения', 'гемолитическая анемия'],
             ['CO[C@@]1(NC(=O)Cc2cccs2)C(=O)N2C(C(=O)O)=C(COC(N)=O)CS[C@@H]21', 'лейкопения', 'гранулоцитопения', 'нейтропения', 'анемия', 'томбоцитопения', 'супрессия костного мозга', 'гемолитическая анемия'],
             ['CCN1CCN(C(=O)N[C@@H](C(=O)N[C@@H]2C(=O)N3C(C(=O)O)=C(CSc4nnnn4C)CS[C@H]23)c2ccc(O)cc2)C(=O)C1=O', 'тромбоцитопения', 'гипопротромбинемия'],
             ['CO/N=C(\C(=O)N[C@@H]1C(=O)N2C(C(=O)[O-])=C(C[n+]3cccc4c3CCC4)CS[C@H]12)c1csc(N)n1', 'аритмия', 'лейкопения', 'гемолитическая анемия', 'тромбоцитопения'],
             ['COCC1=C(C(=O)O)N2C(=O)[C@@H](NC(=O)/C(=N\OC)c3csc(N)n3)[C@H]2SC1', 'тромбоцитоз', 'тромбоцитопения', 'лейкопения', 'эозинофилия'],
             ['CC(C)(O/N=C(\C(=O)N[C@@H]1C(=O)N2C(C(=O)[O-])=C(C[n+]3ccccc3)CS[C@H]12)c1csc(N)n1)C(=O)O', 'гемолитическая анемия'],
             ['CCO/N=C(\C(=O)N[C@@H]1C(=O)N2C(C(=O)[O-])=C(Sc3nc(-c4cc[n+](C)cc4)cs3)CS[C@H]12)c1nsc(NP(=O)(O)O)n1', 'анемия', 'лейкопения', 'нейтропения', 'тромбоцитопения', 'эозинофилия'],
             ['CO/N=C(\C(=O)N[C@@H]1C(=O)N2C(C(=O)O)=C(CSc3nc(=O)c(O)nn3C)CS[C@H]12)c1csc(N)n1', 'анемия', 'гемолитическая анемия', 'гранулоцитопения'],
             ['CO/N=C(\C(=O)N[C@@H]1C(=O)N2C(C(=O)O)=C(COC(N)=O)CS[C@H]12)c1ccco1', 'эозинофилия', 'нейтропения', 'лейкопения', 'анемия', 'гемолитическая анемия', 'тромбоцитопения', 'агранулоцитоз', 'гипопротромбинемия'],
             ['C[C@@H](O)[C@H]1C(=O)N2C(C(=O)O)=C(S[C@@H]3CN[C@H](C(=O)N(C)C)C3)[C@H](C)[C@H]12', 'агранулоцитоз'],
             ['C[C@@H](O)[C@H]1C(=O)N2C(C(=O)O)=C(S[C@@H]3CN[C@H](C(=O)Nc4cccc(C(=O)O)c4)C3)[C@H](C)[C@H]12', 'анемия', 'тромбоцитоз', 'тромбофлебит', 'нейтропения'],
             ['C[C@H]1[C@H](NC(=O)/C(=N\OC(C)(C)C(=O)O)c2csc(N)n2)C(=O)N1S(=O)(=O)O', 'панцитопения', 'нейтропения', 'тромбоцитопения', 'анемия', 'эозинофилия', 'лейкоцитоз', 'тромбоцитоз'],
             ["CO/N=C(\C(=O)N[C@@H]1C(=O)N2C(C(=O)O)=CCS[C@H]12)c1csc(N)n1", "эозинофилия", "тромбоцитоз", "анемия", "гемолитическая анемия", "лейкопения", "нейтропения", "агранулоцитоз", "тромбоцитопения", "панцитопения"],
             ["CO/N=C(\C(=O)N[C@@H]1C(=O)N2C(C(=O)[O-])=C(C[N+]3(C)CCCC3)CS[C@H]12)c1csc(N)n1", "анемия"],
             ["C=CC1=C(C(=O)O)N2C(=O)[C@@H](NC(=O)/C(=N\OCC(=O)O)c3csc(N)n3)[C@H]2SC1", "эозинофилия", "лейкопения", "тромбоцитопения", "нейтропения", "гемолитическая анемия"],
             ["CO/N=C(\C(=O)N[C@@H]1C(=O)N2C(C(=O)O)=C(COC(C)=O)CS[C@H]12)c1csc(N)n1", "нейтропения", "лейкопения"],
             ["C[C@@H](O)[C@H]1C(=O)N2C(C(=O)O)=C(S[C@@H]3CN[C@H](CNS(N)(=O)=O)C3)[C@H](C)[C@H]12", "анемия", "тромбоцитопения", "лейкопения", "нейтропения"]
             ]


def get_effect_name(index):
    return side_effects[index]

def getDB():
    DB = []

    for i in range(n):
        array = [to_fingerprint(compounds[i][0])]

        for j in range(m):
            if side_effects[j] in compounds[i]:
                array.append(1)
            else:
                array.append(0)

        DB.append(array)

    return DB


print(getDB())


def getLabels():
    array = list()
    for i in range(n):
        array.append(int(to_fingerprint(compounds[i][0]), 2))
    return array


print(getLabels())


def getSideEffects():
    effects = list()
    for i in range(n):
        array = list()
        for j in range(m):
            if side_effects[j] in compounds[i]:
                array.append(1)
            else:
                array.append(0)

        effects.append(array)
    return effects


print(getSideEffects())


def getTestLabels():
    array = list()
    for i in range(n, len(compounds)):
        array.append(int(to_fingerprint(compounds[i][0]), 2))
    return array


print(getTestLabels())


def getTestSideEffects():
    effects = list()
    for i in range(n, len(compounds)):
        array = list()
        for j in range(m):
            if side_effects[j] in compounds[i]:
                array.append(1)
            else:
                array.append(0)
        effects.append(array)
    return effects


print(getTestSideEffects())
