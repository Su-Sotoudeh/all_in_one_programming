path = '/home/su/Desktop/RepeatMasker2.lib'
path2 = '//media/su/Bioinf/Dfam38/dfam-fasta-download.fasta'
bird_keywords = [
    "Aves", "Galliformes", "Passeriformes", "Anseriformes", "Accipitriformes",
    "Falconiformes", "Columbiformes", "Psittaciformes", "Strigiformes",
    "Apodiformes", "Gruiformes", "Pelecaniformes", "Phoenicopteriformes",
    "Ciconiiformes", "Coraciiformes", "Caprimulgiformes", "Cuculiformes",
    "Charadriiformes", "Procellariiformes", "Sphenisciformes", "Struthioniformes",
    "Tinamous", "Rheiformes", "Casuariiformes", "Apteryx", "Passer", "Columba",
    "Gallus", "Anas", "Falco", "Buteo", "Corvus", "Taeniopygia", "Melopsittacus",
    "Zebra_finch", "Chicken", "Goose", "Duck", "Finch", "Parrot", "Owl", "Penguin"
]
TEs_list = []
counter = 0
from Bio import SeqIO
for TE in SeqIO.parse(path2, "fasta"):
    taxson  =TE.description.split()[1]
    if taxson in TEs_list:
        continue
    TEs_list.append(taxson)
    counter +=1
    #if taxson in bird_keywords:
    print(counter,taxson)