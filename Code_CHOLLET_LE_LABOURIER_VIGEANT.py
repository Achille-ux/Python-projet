##Projet en groupe
#
#Chollet Héloïse Le Labourier Margot Vigeant Valentin

##Imports

import urllib.request
import ssl
import csv
import os

##Variable globale
#Le dictionnaire CODEAA permet grâce à la fonction traduction retrouvée plus tard de convertir les AA depuis le code 1 lettre vers le code 3 lettres
#L'attribut est l'AA en code une lettre et sa clé ets l'AA en code 3 lettres

CODEAA= {"A":"Ala", "R":"Arg", "D":"Asp", "N":"Asn", "C":"Cys", "E":"Glu", "Q":"Gln", "G":"Gly", "H":"His", "I":"Ile", "L":"Leu", "K":"Lys", "M":"Met", "F":"Phe", "P":"Pro", "S":"Ser", "T":"Thr", "W":"Trp", "Y":"Tyr", "V":"Val"}

##Fonctions :

"""Une fonction qui permet de donner le nombre de numeros d'accession"""

def NumAccession():
    listeAcces=[] #Création d'ue liste vide dans laquelle seront stockés les numéros d'accession des protéines à analyser
    num=input("Avez-vous plus d'un numéro d'accession ? ")
    while num.upper()!="OUI" and num.upper()!="NON": #La méthode upper() renvoie la chaîne en majuscules. Ici l'entrée dans la boucle while ne se fait que si la réponse de l'utilisateur est oui ou non
        print("Desole je ne comprend pas la reponse. Veuillez repondre par oui ou non. ") #Si l'utilisateur a entré un autre terme que oui ou non
        num=input("Avez-vous plus d'un numéro d'accession ? ") #La question est reposée
    if num.upper()=="OUI": #Si l'utilisateur a plusieurs numéros d'accession
        fichier=input("entrez le chemin d'accès du fichier: ")#On récupère les numéros d'accession stockés dans un fichier
        try:
            fich=open(fichier,"r") #La fonction open ouvre le fichier et le renvoie en mode lecture
            li=fich.readline() #la fonction .readline lit une ligne entière du fichier
            while li !="": #Tant que les lignes du fichier ne sont pas vides alors on parcourt le fichier
                if li[-1:]=="\n": #Quand le programme a fini de lire une ligne, marquée par un retour chariot
                    listeAcces.append(li[:-1]) #Alors on ajoute à la liste le numéro d'accession lu à la ligne précédente
                else:
                    listeAcces.append(li) #Le dernier numéro d'accession a été ajouté à la liste
                li=fich.readline()
        except:
            print("Le chemin d'accès est invalide") #L'utilisateur a fait une erreur dans le chemin d'accès, le fichier n'est pas trouvable par ce chemin.
            fichier=input("entrez le chemin d'accès du fichier: ") #On lui repose la question
    elif num.upper()=="NON" : #Si l'utilisateur n'a qu'un seul numéro d'accession
        fi=input("entrez le numéro d'accession de la fiche: ")
        listeAcces.append(fi) #On ajoute à la liste jusqu'à présent vide nommée listeAcces le numéro d'accession de la protéine
    return (listeAcces, num)

""" Une fonction qui permet de recuperer dans une liste les informations de la fiche Uniprot trouvee sur Internet"""

def recupfiche(listeAcces):
    l=[] #Création d'une liste vide
    ls=[] #Les élèments de la liste vide ls seront ensuite stockés dans la liste l
    for i in range (0, len(listeAcces)): #Pour i compris entre le premier et le dernier élèment de la liste
        try:
            context = ssl._create_unverified_context() #ssl._create_unverified_context permet de faire une vérification des certificats présentés par les serveurs HTTPS
            u=urllib.request.urlopen("https://rest.uniprot.org/uniprotkb/"+listeAcces[i].upper()+".txt", context=context) #urlopen permet d'ouvrir le lien auquel on a associé le numéro de la protéine pour arriver sur sa fiche PDB
            uniprotlines=u.readlines()#à uniprotlines est associé l'ensemble des lignes lues dans la fiche UniProt
            u.close()
        except: #S'il y a une erreur dans l'url
            print("Probleme lors de la lecture du fichier")
            print("https://rest.uniprot.org/uniprotkb/"+listeAcces[i].upper()+".txt") #Affichage du lien généré grâce aux informations entrées par l'utilisateur pour qu'il visualise sont erreur
        else:
            for ligne in uniprotlines:
                ls.append (ligne.decode("utf8").strip()) #Un décodage du contenu web est nécessaire (UTF8 : Unicode Transformation Format avec 8 bits utilisés dans l’endocodage).
            l.append(ls) #Les données stockées dans la sous-liste ls sont ajoutées à la liste l
        i=i+1
        ls=[]
    return (l)

""" Une fonction qui permet de creer un dictionnaire dont les elements sont les elements de la liste que l'on souhaite recuperer"""

def triinformations(l):
    nomproteine, basedonnes, identifiant, longueur, gene, poids, description2, espece2, sequenceproteique, codeP, codeI, codePro = "", "", "", "", "", "", "", "", "", "", "", ""
    p=0
    listeDicos=[] # C'est dans cette liste finale que nous allons stocker toutes les informations pertinentes sur la protéine à rendre à l'utilisateur
    for el in l: #Pour chaque élèment de la liste l
        dicoelements={}
        description, espece, sequence=[], [], [],
        listeDicos.append(el) #A la liste est ajouté l'élément el
        n=0   #On initialise un compteur pour compter les domaines transmembranaires
        i=0
        while i<len(el): #Tant que l'on n'a pas parcouru tout l'élèment el
            e=el[i].split() #La fonction split() permet de séparer chaque élèment de la ligne comme un élèment unique à partir de l'élèment ligne globale initiale
            if e[0]=="ID": #L'acronyme ID est le premier élèment de la ligne
                nomproteine=nomproteine+e[1] #A nom de la protéine est associé l'élèment 1 de la ligne ID qui correspond au nom de la protéine
                longueur=longueur+e[3] #A longueur est associé l'élèment 3 de la ligne
                basedonnes=basedonnes+(e[2])
                basedonnes=basedonnes.replace(";", "")#On supprime les points virgules
            elif e[0]=="AC": #IDEM que popur ID, l'accronyme AC est désigné comme premier élèment de la liste
                identifiant=identifiant+e[1]
                identifiant=identifiant.replace(";", "")
            elif e[0]=="GN":
                h=e[1].split("=") #On sépare les caractères séparés par un =
                gene=gene+h[1]
                gene=gene.replace(";", "")
            elif e[0]=="OS":
                h=el[i].split("  ")
                espece=espece+h[1:]
                espece2=' '.join(espece)#On joint les différentes lignes avec OS
            elif e[0]=="SQ":
                h=el[i].split("  ")
                poids=poids+h[3]
                poids=poids.replace("MW;", "") #La fonction replace permet ici de supprimer "Molecular Weight"
            elif len(e[0])>2: #Lorsuqe que la longueur de l'élèment est supérieure à 2, cela signifie qu'il s'agit de la séquence protéique selon la configuration habituelle d'une fiche PDB
                sequence=sequence+e
                sequenceproteique=''.join(sequence)#La fonction join permet ici d'ajouter un espace entre les élèments
                sequenceproteique=sequenceproteique.replace(" ","")
            elif e[0]=="DE":
                description=description+e[1:]
                description2=' '.join(description)
                description2=description2.replace(";", "")
            elif e[0]=="DR": #Dans l'élèment DR, on trouve toutes les références PDB, Interpro et Prosite
                if e[1]=="PDB;":
                    codeP=codeP+e[2]
                    codeP=codeP.replace(";", "-") #On remplace les points virgules par des tirets
                elif e[1]=="InterPro;":
                    codeI=codeI+e[2]
                    codeI=codeI.replace(";", "-")
                elif e[1]=="PROSITE;":
                    codePro=codePro+e[2]
                    codePro=codePro.replace(";", "-")
            elif e[0]=="FT" and e[1]=="TRANSMEM": #Lorsque l'acronyme FT et TRANSMEM sont reconnu comme élèments 1 et 2 sur la ligne cela fait référence aux domaines transmembranaires
                n=n+1 #Le compteur est incrémenté de 1 à chaque domaine transmembranaire dénommé par "FT" et "TRANSMEM"
            i=i+1
        dicoelements={"Nom de la proteine":nomproteine, "Base de donnes" : basedonnes, "Identifiant": identifiant, "Longueur de la sequence (AA)": longueur, "Nom du gene": gene, "Espece": espece2, "Poids moleculaire": poids, "Description": description2, "References PDB":codeP[:-1], "References InterPro": codeI[:-1], "References PROSITE":codePro[:-1], "Nombre de domaines transmembranaires":n, "Sequence proteique": sequenceproteique}  #On ajoute :-1 pour enlever le tiret de la fin
        #Toutes les données obtenues dans la boucle while précédente sont stockées dans le dictionnaire dicoelement
        if dicoelements["References PDB"]=="": #S'il n'y a pas de référence PDB dans la fiche
            dicoelements["References PDB"]="Pas de donnees"
        if dicoelements["References InterPro"]=="": #S'il n'y a pas de référence InterPro dans la fiche
            dicoelements["References InterPro"]="Pas de donnees"
        if dicoelements["References PROSITE"]=="": #S'il n'y a pas de référence Prosite dans la fiche
            dicoelements["References PROSITE"]="Pas de donnees"
        if dicoelements["Nom du gene"]=="":
            dicoelements["Nom du gene"]="Pas de donnes"
        listeDicos[p]=dicoelements
        p=p+1 #On passe à l'élèment suivant de la liste l
        dicoelements={} #le dictionnaire est vidé pour être rempli et analysé par l'élèment suivant de la liste l
        nomproteine, basedonnes, identifiant, longueur, gene, poids, description2, espece2, sequenceproteique, codeP, codeI, codePro = "", "", "", "", "", "", "", "", "", "", "", ""
    return (listeDicos)

""" Une fonction qui permet de mettre toutes ces informations sous la forme d'un tableau"""

def tableau(listeDicos):
    nom=input("Comment souhaitez-vous appeler le fichier ? ")
    with open(nom+".csv", "w") as csvfile: #Csvfile permet de transformer un fichier texte brut en organisant les données sous forme de tableau
        writer=csv.DictWriter(csvfile, fieldnames=["Nom de la proteine", "Base de donnes", "Identifiant", "Longueur de la sequence (AA)", "Nom du gene", "Espece", "Poids moleculaire", "Description", "References PDB", "References InterPro", "References PROSITE", "Nombre de domaines transmembranaires", "Sequence proteique"], delimiter=";")
        #DictWriter est une fonction qui crée un objet qui opère comme un transcripteur ordinaire mais qui produit les lignes de sortie depuis des dictionnaires
        writer.writeheader() # Write.header est une fonction qui écrit une ligne contenant les noms de champs
        for i in range (0, len(listeDicos)):
            writer.writerow(listeDicos[i]) #write.row va écrire le paramètre row vers le fichier associé au transcripteur, une ligne pour une protéine
    csvfile.close()

"""Une fonction qui permet d'ecrire dans un fichier la sequence multifasta"""

def ecrirefichierfasta(listeDicos):
    i=0
    fichierfasta=open("fichier_multifasta.txt", "w") #La fonction open permet d'ouvrir la ficher suivant dans la parenthèse en mode écriture
    for d in range (0, len(listeDicos)): #Dans cette boucle, chaque élèment de la listeDico est parcouru
        fichierfasta.write("> {} {} {} {} {}".format(listeDicos[d]["Nom de la proteine"], listeDicos[d]["Identifiant"], listeDicos[d]["Base de donnes"], listeDicos[d]["Longueur de la sequence (AA)"], listeDicos[d]["Nom du gene"]))
        while i<len(listeDicos[d]["Sequence proteique"]):
            if i<len(listeDicos[d]["Sequence proteique"]):
                fichierfasta.write("\n {}".format(listeDicos[d]["Sequence proteique"][i:i+80])) #On écrit 80 caractères par ligne
                i=i+80
            else:
                fichierfasta.write("\n {}".format(listeDicos[d]["Sequence proteique"][i:])) #Si il reste moins de 80 caractères, on écrit tout ce qui reste
        fichierfasta.write("\n")
        i=0
    fichierfasta.close()

"""Une fonction qui traduit en code a 3 lettres"""

def traduction (listeDicos):
    seq3lettres=[] #Liste vide dans laquelle sera stockée la séquence à trois lettres en acides aminés de la protéine
    for d in range (0, len(listeDicos)):
        s3l=""
        i=0
        while i<len(listeDicos[d]["Sequence proteique"]): #ON parcourt la sequence
            s3l=s3l+CODEAA[listeDicos[d]["Sequence proteique"][i]] #A s3l est attribué le nom de l'acide aminé en code 3 lettres stocké dans la listeDico et correspondant à l'élèment sequence proteique
            i=i+1
        seq3lettres.append(s3l) #Un acide aminé est ajouté par tour à la liste seq3lettres
        s3l="" #On remet à vide la valeur de s3l pour y traduire le prochain acide aminé de la chaine protéique
    return (seq3lettres)

"""Une fonction qui compte chaque acide aminé"""

def analysesequence(seq3lettres, listeDicos):
    listeCompteur=[] #Liste vide dans laquelle seront stokées les occurences de chaque acide aminé
    listeFreq=[] #Liste vide dans laquelle seront stokées les fréquences de chaque acide aminé
    for d in range(0, len(seq3lettres)):
        compteur={"Ala":0, "Arg":0, "Asp":0, "Asn":0, "Cys":0, "Glu":0, "Gln":0, "Gly":0, "His":0, "Ile":0, "Leu":0, "Lys":0, "Met":0, "Phe":0, "Pro":0, "Ser":0, "Thr":0, "Trp":0, "Tyr":0, "Val":0}
        frequence={"Ala":0, "Arg":0, "Asp":0, "Asn":0, "Cys":0, "Glu":0, "Gln":0, "Gly":0, "His":0, "Ile":0, "Leu":0, "Lys":0, "Met":0, "Phe":0, "Pro":0, "Ser":0, "Thr":0, "Trp":0, "Tyr":0, "Val":0}
        i=0
        while i<len(seq3lettres[d]):#Boucle qui parcours toute la séquence en AA 3 lettres
            compteur[seq3lettres[d][i:i+3]]=compteur[seq3lettres[d][i:i+3]]+1 #A chaque AA rencontré, son compteur est incrémenté de 1
            frequence[seq3lettres[d][i:i+3]]=frequence[seq3lettres[d][i:i+3]]+1 #A chaque AA rencontré, son compteur est incrémenté de 1
            i=i+3
        for j in frequence:
            frequence[j]=((frequence[j])/int(listeDicos[d]["Longueur de la sequence (AA)"]))*100 #Permet le calcul de la fréquence en prenant l'occurence de chaque AA divisé par la longueur de la séquence totale. Puis multiplier par 100 pour l'avoir en pourcentage
        sorted_items = sorted(compteur.items(), key=lambda x: -x[1])  #Pour que ça trie dans l'odre croissant
        compteur.clear()
        compteur.update(sorted_items)
        listeCompteur.append(compteur) #L'AA rencontré est additionné à chaque tour de boucle à la liste compteur
        listeFreq.append(frequence)
    return (listeCompteur, listeFreq)

"""Une fonction qui affiche le fichier contenant le nom de l'acide amine en code 3 lettres, le nombre de fois qu'il apparait et sa frequence"""

def fichieranalysesequence(listeCompteur, listeDicos, listeAcces, listeFreq):
    for d in range(0, len(listeCompteur)):
        fichieranalyse=open("analyse_sequence"+listeAcces[d]+".txt", "w") #Création d'un fichier au format txt en mode écriture
        fichieranalyse.write("AC={}".format(listeDicos[d]["Identifiant"])) #On écrit l'identifiant de la protéine dans le fichier
        for i in listeCompteur[d]:
            fichieranalyse.write("\n{} {} {:.2f}".format(i, listeCompteur[d][i], listeFreq[d][i]))#Ajout de la fréquence et de l'occurence en AA de la protéine
        fichieranalyse.close()

##Programme principal

print("Recuperation d'une ou plusieurs fiche(s) UNIPROT ")
listeAcces, num=NumAccession()
l=recupfiche(listeAcces)
listeDicos=triinformations(l)
if num.upper()=="NON": #Dans le cas ou on n'a qu'un seul numéro d'accession
    print("Nom de la proteine recherchee : {}\nBase de donnees : {}\nIdentifiant de la proteine : {}\n\nLongueur de la sequence : {} acides amines\nNom du gene : {}\nEspece : {}\nPoids moléculaire de la proteine : {}\n\nDescription : {}\n\nReferences PDB : {}\nReferences PROSITE : {}\nReferences InterPro : {}\n\nNombre de domaines transmembranaires : {}\n\nSequence Proteique : {}.".format(listeDicos[0]["Nom de la proteine"], listeDicos[0]["Base de donnes"], listeDicos[0]["Identifiant"], listeDicos[0]["Longueur de la sequence (AA)"], listeDicos[0]["Nom du gene"], listeDicos[0]["Espece"], listeDicos[0]["Poids moleculaire"], listeDicos[0]["Description"], listeDicos[0]["References PDB"], listeDicos[0]["References PROSITE"], listeDicos[0]["References InterPro"], listeDicos[0]["Nombre de domaines transmembranaires"], listeDicos[0]["Sequence proteique"] )) #On affiche dans la console les informations
    qu=input("Voici les informations. Voulez-vous les enregistrer dans un fichier ? ")
    while qu.upper()!="OUI" and qu.upper()!="NON": #Tant que l'utilisateur ne répond pas la réponse attendue soit oui ou non
        print("Desole je ne comprend pas la reponse, veuillez repondre par oui ou non")
        qu=input("Voici les informations. Voulez-vous les enregistrer dans un fichier ? ") #On lui repose la question
    if qu.upper()=="OUI":
        tableaufinal=tableau(listeDicos)
    print("\n> {} {} {} {} {}\n{}".format(listeDicos[0]["Nom de la proteine"], listeDicos[0]["Identifiant"], listeDicos[0]["Base de donnes"], listeDicos[0]["Longueur de la sequence (AA)"], listeDicos[0]["Nom du gene"], listeDicos[0]["Sequence proteique"]))
    qi=input("\nVoici la sequence. Voulez vous l'enregistrer dans un fichier au format fasta ? ")
    while qi.upper()!="OUI" and qi.upper()!="NON":
        print("Desole je ne comprend pas la reponse, veuillez repondre par oui ou non")
        qi=input("\nVoici la sequence. Voulez vous l'enregistrer dans un fichier au format fasta ? ")
    if qi.upper()=="OUI":
        fichierfasta=ecrirefichierfasta(listeDicos)
        print("Fichier enregistre")
    qr=input("\nVoulez-vous une analyse de la sequence d'acides amines ? ")
    while qr.upper()!="NON" and qr.upper()!="OUI":
        print("Desole je ne comprend pas la reponse, veuillez repondre par oui ou non")
        qr=input("Voulez-vous une analyse de la sequence ? ")
    if qr.upper()=="OUI":
        seq3lettres=traduction(listeDicos)
        listeCompteur, listeFreq=analysesequence(seq3lettres, listeDicos)
        fichiersequence=fichieranalysesequence(listeCompteur, listeDicos, listeAcces, listeFreq)
else: #Dans le cas où on a plusieurs numéros d'accession
    tableaufinal=tableau(listeDicos)
    print("Les informations des fiches sont enregistrees dans un tableau. ")
    qi=input("Voulez vous enregistrer les sequences dans un fichier au format fasta ? ")
    while qi.upper()!="OUI" and qi.upper()!="NON":
        print("Desole je ne comprend pas la reponse, veuillez repondre par oui ou non")
        qi=input("\nVoulez vous enregistrer les sequences dans un fichier au format fasta ? ")
    if qi.upper()=="OUI":
        fichierfasta=ecrirefichierfasta(listeDicos)
        print ("Fichier enregistré")
    qr=input("\nVoulez-vous une analyse de la sequence d'acides amines ? ")
    while qr.upper()!="NON" and qr.upper()!="OUI":
        print("Desole je ne comprend pas la reponse, veuillez repondre par oui ou non")
        qr=input("Voulez-vous une analyse de la sequence ? ")
    if qr.upper()=="OUI":
        seq3lettres=traduction(listeDicos)
        listeCompteur, listeFreq=analysesequence(seq3lettres, listeDicos)
        fichiersequence=fichieranalysesequence(listeCompteur, listeDicos, listeAcces, listeFreq)
input("Vos fichiers sont dans le dossier ou se trouve ce programme. Tapez sur entree pour terminer")
