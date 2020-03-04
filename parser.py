import re # import per l'utilizzo delle regex
import sys # import per chiudere il programma
from urllib.parse import urlparse

"""
    TODO: compressed file
    
    righe che sono univoche: fileFormat, 'fileDate', 'reference'

"""




# Definisco i pattern per le regex comuni 
id_pattern = "ID=[a-zA-Z0-9_]+"
desc_pattern = "Description=\"[^\"]+\""
num_patter = "Number=(?P<num>([0-9]+|A|R|G|.))"




"""
    funzione che legge tutto il contenuto di un file   
    filename = nome del file
    compressed = specifica se il file è in formato compresso oppure no
"""
def reader(filename, compressed = False):
    # TODO: gestione della compressione
    
    # apro e leggo tutto il contenuto del file
    with open(filename, 'r') as vcf_file:
        lines = vcf_file.readlines()

    vcf_file.close() # chiudo il file
    return lines # ritorno il contenuto del file
    
    
    
"""
   funzione che controlla che la prima riga sia del tipo 
   ##fileformat=VCFvNumber.Number
   Return 0 = va bene
   -1 = errore nel formato
"""
def fileformat(line):
    response =  re.search("##fileformat=VCFv[0-9].[0-9]" ,line, re.IGNORECASE) # controllo che la prima riga del file vada bene
    if (response):
        return 0
    else:
        return -1
    
   
"""
   controllo la riga FILTER
   ##filter=<ID=ID,Description="description">
   Return 0 = va bene
   -1 = errore nel formato
""" 
def myFilter(line):
    response = re.search("FILTER=<" + id_pattern + ",\s*" + desc_pattern + ">" ,line, re.IGNORECASE)
    if (response):
        return 0
    else:
        return -1


"""
   controllo la riga info
   ##info=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">
   Return 0 = va bene
   -1 = errore nel formato
   2 = warning, numero diverso da 0 con type = flag!
""" 
def myInfo(line):
    response = re.search("INFO=<" + id_pattern + ",\s*" + num_patter + ",\s*" +
                     "Type=(?P<type>(Integer|Float|Flag|Character|String)),\s*" + 
                     desc_pattern + "(,\s*Source=(\")*[a-zA-Z0-9_\s]+(\")*)*" +
                     "(,\s*Version=(\")*[0-9]+(\")*)*>" ,line, re.IGNORECASE)
    
    # se ho matchato prima di ritornare faccio il controllo incrociato Flag-Number
    if (response):
        # controllo che se type = flag allora number = 0.. E' UNO SHOULD BE, QUINDI NON OBBLIGATORIO
        if (response.group("type") == "Flag" and response.group("num") != "0"):
            return 2
        return 0
    else:
        return -1


"""
   controllo la riga format
   ##info=<ID=ID,Number=number,Type=type,Description="description">
   Return 0 = va bene
   -1 = errore nel formato
""" 
def myFormat(line):
    response = re.search("FORMAT=<" + id_pattern + ",\s*" + num_patter + ",\s*" +
                     "Type=(?P<type>(Integer|Float|Character|String)),\s*" + 
                     desc_pattern + ">", line, re.IGNORECASE)
    
    # se ho matchato prima di ritornare faccio il controllo incrociato Flag-Number
    if (response):
        return 0
    else:
        return -1


"""
   controllo la riga alt
   ##info=<ID=ID,Description="description">
   Return 0 = va bene
   -1 = errore nel formato
""" 
def myAlt(line):
    
    response = re.search("ALT=<" + id_pattern + ",\s*" + desc_pattern + ">", line, re.IGNORECASE)
    
    # se ho matchato prima di ritornare faccio il controllo incrociato Flag-Number
    if (response):
        return 0
    else:
        return -1

    
    
    
    
    
    
"""
    Funzione che data una stringa ripulita dagli ## iniziali controlla che 
    sia well-formed e nel caso di righe che devono essere univoche controlla
    che sia la prima volta che la si incontra!
"""
def metadataParser(line):
    """
        1) cerco il primo = dato che ogni stringa è una coppia key-value
        2) faccio uno "switch" per capire in quale metadato mi trovo
    """
    indexEqual = line.find("=") 
    
    info = line[:indexEqual]
    info = info.upper() # metto la riga in maiuscolo, cosi matcho anche con le righe in minuscolo
    result = -1
    if (info == "FILTER"):
        result = myFilter(line)
    elif (info == "INFO"):
        result = myInfo(line)
    elif (info == "FORMAT"):
        result = myFormat(line)
    elif (info == "ALT"):
        result = myAlt(line)
    else:
        result = 0
    return result



# Nel parser controllo subito che la prima riga sia il fileFormat, se lo è continuo
def parser(allLines):
    # controllo la prima riga del file!
    result = fileformat(allLines[0]) 
    if (result == -1): # se da errore chiudo il programma e stampo un messaggio
        sys.exit("SyntaxError at line: " + str(0))
    else:
        #del allLines[0] # cancello la prima riga, l'ho già parsata
        
        #totalLines = len(allLines) # salvo il numero totale di righe lette
        actualLine = 2 # riga attuale di lettura
        metadataAllows = True # variabile booleana che mi dice se sono ancora nella sezione ## oppure sono nella sezione dati
        for line in allLines[1:] :
            
            isMetadata = re.search("##", line) ## controllo se è una riga dati o di metadati
            
            # Se ho trovato una stringa di metadati e sono nella sezione metadati allora:
            # tolgo i ## iniziali e passo la stringa al parser ad hoc per i metadati
            if (isMetadata and metadataAllows):
                result = metadataParser(line[2:])
                
            #return could be 0: OK, -1: error, 2: special warning for flag type
            
            if(result == -1):
                sys.exit("SyntaxError at line: " + str(actualLine))
                
            elif(result == 2): #special warning for Flag type in format line
                print("Warning: number should be zero with Flag type at line: " + str(actualLine))
                            
            actualLine = actualLine + 1 
            # se la riga che ho letto va bene aumento di uno il contatore delle righe
            #actualLine = actualLine + 1 
    

    
    

filename = './Test/example-4.2.vcf'

lines = reader(filename)
parser(lines)