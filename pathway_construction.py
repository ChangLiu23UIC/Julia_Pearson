import pandas as pd

# Your lists of elements
AURKA = ['Cdc25', 'Aurora-A', 'PP2A', 'Xkid', 'G', 'CPEB', 'Calm', 'PKA', 'GPCR', 'CycB4', 'Emi1', 'Cdc20', 'CycB',
         'Rsk1', 'C3H-4', 'Rec8', 'Ca2+', 'Cdk2', '14-3-3', 'SCF', 'B56', 'Plk1', 'PKB', 'Separase', 'Mad1', 'APC/C',
         'Hsp90', 'MAPK12', 'Mad3', 'Akt', 'INS', 'p42MAPK', 'MAPK', 'CaMKII', 'IGF-1', 'Raf', 'Mos', 'Eg2',
         'Plkk1', 'Seo', 'CycA', 'IP3R', 'Sgo', 'PIP2', 'SMC1', 'CaN', 'PI3K', 'PGR',
         'Bub1', 'XPR-1', 'Mad2', 'AC', 'Plx1', 'Ras', 'STAG3', 'IGF1R', 'β-TRCP', 'Insulin', 'Cdc2', 'Securin',
         'Cdc25C', 'SMC3', 'Fizzy', 'Emi2', 'MPF', 'IP3', 'CycE', 'PP1', 'Ringo', 'PDE3', 'Rsk2', 'Myt1', 'PIP3',
         'AR', 'CycB1', 'MEK1']

AURKB = ['RBBP4', 'TRIP13', 'Smad2', 'FBRS', 'SKP1', 'PHC1', 'SCML1', 'Cip1', 'TEX10', 'E2F1,2,3', 'PHC3',
               'ASXL', 'Rad21', 'Smad3', 'CBX8', 'SCF', 'Mad1', 'Mdm2', 'SCMH1', 'Nde80', 'PP2A', 'CBX4', 'SUZ12',
               'NIPBL', 'RNF2', 'Cdh1', 'E2F6', 'ATM', 'ATR', 'p53', 'UBE2D', 'USP16', 'MAX', 'PHF1', 'PTTG', 'HDAC8',
               'SFMBT', 'AEBP2', 'MCM', 'BubR1', 'MBD5', 'Securin', 'Cdt1', 'Cdc45', 'CDK7', 'Cdc6', 'HDAC', 'PDS5',
               'Separin', 'APC/C', 'PCGF5', 'CDK2', 'Sororin', 'STAG2', 'Ink4c', 'AUTS2', 'STAG1', 'LCOR', 'Bub3',
               'Cdc7', 'Mps1', 'YAF2', 'CDK4,6', 'Dbf4', 'RYBP', 'EZH2', 'PCGF3', 'TFDP1', 'p57', 'Cdc25B,C', 'FOXK',
               'CBX2', 'PCGF6', 'Ink4d', 'p300', 'E2F5', 'CBX6', 'DNA-PK', 'BAP1', 'Cdc20', 'EFD', 'E2F4', 'GADD45',
               'TGFβ', 'Kip1,2', 'YY1', 'RBBP7', 'DP-1', 'MTBP', 'SCML2', 'AURK', 'ARF', 'CMT2', 'Sme3', 'ESCO2',
               'Sme1', 'WDR5', 'ATRX', 'ORC', 'HCFC1', 'p107', 'Stag1', 'Skp2', 'AuroraB', 'PCGF4', 'Wee1', 'PCNA',
               'GSK3β', 'Ink4a', 'CycA', 'PCGF2', 'Mad2', 'Emi1', 'PCGF1', 'DCAF7', 'Cdc25A', 'CBX7', 'CSNK2E',
               'Smad4', 'CBX3', 'MGA', 'CDK1', 'EPOP', 'L3MBTL', 'Sgo1', 'EZH1', 'MBD6', 'DDX11', 'PHF19',
               'CycH', 'KDM2', 'RNF1', 'Plk1', 'p21', 'CycD', 'PHC2', 'SMC3', 'SMC1', 'Esp1', 'c-Myc', 'OGT',
               'ESCO1', 'p27', 'Mcl-1', 'BCOR', 'USP7', 'KNL1', 'CycE', 'TICRR', 'Stag2', 'MAU2', 'Bub1', 'MTF2',
               'Chk1,2', 'Wapl', 'JARID2', 'DP-2', 'Miz1', 'Ink4b']

# Make sure both lists have the same length
max_length = max(len(AURKA), len(AURKB))

# Pad the shorter list with None
AURKA.extend([None] * (max_length - len(AURKA)))
AURKB.extend([None] * (max_length - len(AURKB)))

# Create a DataFrame
df = pd.DataFrame({
    "AURKA": AURKA,
    "AURKB": AURKB
})

# Save the DataFrame to an Excel file
df.to_excel("Protein.xlsx", index=False)

# Display the DataFrame
print(df)
