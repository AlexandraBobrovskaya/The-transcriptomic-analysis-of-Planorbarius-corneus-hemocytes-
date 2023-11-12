"""
Suggested usage: 
Feed blast XML file and file with hierarchial key-words for search.

Script should extract info on hits with the description containing specified key-words. Probably extract sequences in other files, gather stats and etc.


in-out: streams and files



test XML file path = "~/Desktop/Herzen/Transcriptome/immunityRelatedGenesSearch/preprocBBDuk27/BlastXML_test.xml"
path to full file~/Desktop/Herzen/Transcriptome/immunityRelatedGenesSearch/preprocBBDuk27/preprocBBduk27aoPredProt_BLASTp_NR_NEW.xml

protein TSV path = "~/Desktop/Herzen/Transcriptome/immunityRelatedGenesSearch/preprocBBDuk27/listOfImmunityRelatedGenes.tsv"

test tsv path = "~/Desktop/Herzen/Transcriptome/immunityRelatedGenesSearch/preprocBBDuk27/list_Test.tsv"

"""
from Bio.Blast import NCBIXML
import argparse
import pandas as pd
import csv

import re


# =============================================================================
# Classes
# =============================================================================


class proteinRecord:
    """
    Object storing all the information about found immunity-related protein coded by transcriptome, includes all HSPs for that polypeptide.
    INTENDED USAGE:
        Create object of this class for every protein hit found in XML file.
    """
    
    def __init__(self,
                 queryName,
                 hit):
        self.name = queryName
        self.hits = [hit]
        self.protein = "notSpec"
        self.CDS = "notSpec"
        self.motherTranscript = "notSpec"
        #          hitName,
        #          hitSpecieID,
        #          hitSpecieSciName,
        #          HSPlen,
        #          HSPident,
        #          classificationInfo):
        # self.queryName = queryName
        # self.hitName = [hitName]
        # self.hitTaxNum = [hitSpecieID]
        # self.hitSpecieSciName = [hitSpecieSciName]
        # self.identity = HSPident/HSPlen*100
        # self.HSP = "notSpec"
        # self.protein = "notSpec"
        # self.CDS = "notSpec"
        # self.motherTranscript = "notSpec"
        # self.classificationInfo = classificationInfo
        
        
    def __str__(self):
        hitString = ''
        for hit in self.hits:
            hitString += str(hit) + '\n'
        # return f"Query ID is {self.name}\nHit name is {self.hitName}\nFamily name is {self.classificationInfo.name}\n \
# Class name is {self.classificationInfo.className}"
        return (f"Query ID is {self.name}\n{len(self.hits)} hits are found \n\
The hits are:\n{hitString}")

################################################################################

class proteinToSearch:
    """
    Object representing a certain class of immunity-related proteins, stores all information about it. 
    
    INTENDED USAGE:
        Filled while reading .tsv file with protein to search information.
        Use "searchKeywords" attribite during BlastHit description parsing. After getting a hit mention reference particular 
        proteinToSearch object in proteinRecord object (i.e. this found protein is of this class).
    
    """
    def __init__(self,
                 shortName,
                 fullName,
                 aliases,
                 subclassName,
                 className):
        self.name = shortName
        self.fullName = fullName
        if not aliases:
            self.searchKeywords = [self.name.lower()]
        else:
            self.searchKeywords = [self.name.lower()] + [i.strip() for i in aliases.lower().split(",")]
        # if not aliases:
        #     self.searchKeywords = [self.name]
        # else:
        #     self.searchKeywords = [self.name] + aliases.replace(" ", "").split(",")
        self.subclass = subclassName
        self.className = className
        self.occurrence = 0

    def incrementCounter(self):
        self.occurrence += 1


    def __str__(self):
        """
        

        Returns
        -------
        str
            DESCRIPTION.

        """
        return f'General query record name is {self.name}, it belongs to "{self.className}" class and "{self.subclass}" subclass \n \
            the keyWords for search are {self.searchKeywords}'
            
################################################################################

class BlastHit:
    def __init__(self,
                 hitName,
                 keyword,
                 eValue,
                 hitAccession,
                 hitSpecieID,
                 hitSpecieSciName,
                 HSPlen,
                 HSPident,
                 classificationInfo):
        '''
        Class for storage of a single Blast protein hit.

        Parameters
        ----------
        hitName : TYPE
            DESCRIPTION.
        hitAccession : TYPE
            DESCRIPTION.
        hitSpecieID : TYPE
            DESCRIPTION.
        hitSpecieSciName : TYPE
            DESCRIPTION.
        HSPlen : TYPE
            DESCRIPTION.
        HSPident : TYPE
            DESCRIPTION.
        classificationInfo : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        self.hitName = hitName
        self.keyword = keyword
        self.eValue = eValue
        self.accession = hitAccession
        self.hitTaxNum = hitSpecieID
        self.hitSpecieSciName = hitSpecieSciName
        self.identity = HSPident/HSPlen*100
        self.HSP = "notSpec"
        self.classificationInfo = classificationInfo

    def __str__(self):
        return (
            f"Hit name is {self.hitName}\n\
                It's {self.hitSpecieSciName}\n\
                    Family name is {self.classificationInfo.name}\n\
                        Class name is {self.classificationInfo.className}"
            
            )

################################################################################

class IPShit:
   def __init__(self,
                hitName,
                keyword,
                eValue,
                classificationInfo
                ):
       '''
       Class for storage of a single IPS hit. Now work only with tsv IPS files.
       Though XML contains much more information.

       Parameters
       ----------
       hitName : TYPE
           DESCRIPTION.
       keyword : TYPE
           DESCRIPTION.
       eValue : TYPE
           DESCRIPTION.
       classificationInfo : TYPE
           DESCRIPTION.

       Returns
       -------
       None.

       '''
       self.hitName = hitName
       self.keyword = keyword
       self.eValue = eValue
       self.classificationInfo = classificationInfo
       
   def __str__(self):
    return (
        f"Hit name is {self.hitName}\n\
Family name is {self.classificationInfo.name}\n\
Class name is {self.classificationInfo.className}"           
           )





# =============================================================================
# Functions
# =============================================================================

def readProteinTSV(filepath, *args, **kwargs):
    """
    Function to convert input tsv-file into list to feed it downstream.

    Parameters
    ----------
    filepath : open file ('r')
        tsv-file containing all info about proteins/domains to find.
        Columns:
            Full name	Primary name	Alias	Subclass	Class
    *args : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    searchList : list
        List where each element corresponds to one row of input tsv file.

    """
    searchList = []

    csvReader = csv.reader(filepath, delimiter="\t")
    next(csvReader)                 ###To skip header which is the 1st line of a file

    for row in csvReader:
        searchList.append(proteinToSearch(row[1], row[0], row[2], row[3], row[4]))
        
    # for entry in searchList:
    #     print(entry)
    
    return searchList


################################################################################


def searchBlastXML(XML, queryList, level = "subclass", reportFrequency = 5000):
    """
    Main function checking presence of keywords of interest in corresponding fields of Blast XML file.
    On hit extracts additional info.

    Parameters
    ----------
    XML : XML 2.0 produced by NCBI Blast
        XML file with the results of Blast search.
    queryList : list
                Output of readProteinTSV function
        List describing proteins/domains to look for.

    Returns
    -------
    None.

    """
    
    relevantLevels = {"domain", "subclass", "class"}
    
    if level not in relevantLevels:
        raise ValueError ('Level must be one of the following ("domain", "subclass", "class")')
    elif level == "domain":
        level = "name"
    elif level == "class":
        level = "className"
    
    
    # print(queryList)
    
    counterBlastEntries = 0
    matchesFound = 0
    dictOfSatisfyingProteins = {}
    
    specieNameRegex=re.compile('\[.*\]')    # construction to detect specie name which usually is in the end of a line: [Homo sapiens]
    accessionRegex=re.compile('\|.*\|')    # construction to detect accession number which usually is in the start of a line: |gbp84660.1|

    for entry in NCBIXML.parse(XML):
        matchFound = False
        counterBlastEntries += 1
        # print(f"Polypetide ID is {entry.query}")
        
        numberOfHit = 0 ### Need to recover taxonomy ID cos it's in entry.descriptions.items
        for hit in entry.alignments:
            
# =============================================================================
#             Looking for a specie name in the title because it confuses the script sometimes
# =============================================================================
            
            hitTitle = hit.title.lower()
            
            ### Using regex
            searchSpecieName = specieNameRegex.search(hitTitle)  
            searchAccession = accessionRegex.search(hitTitle)
            if searchSpecieName is not None:
                titleSpecieName = hitTitle[searchSpecieName.start():searchSpecieName.end()]
                hitTitleToSearch = hitTitle[searchAccession.end():searchSpecieName.start()]
                # print(f'Searching in {hitTitleToSearch}')
                # print(titleSpecieName)
            else:
                hitTitleToSearch = hitTitle
            ###
                
            ### Other way of omitting specie name from the searh
            # if '[' in hitTitle:
            #     hitTitleToSearch, titleSpecieName = hitTitle.split('[')
            #     print(titleSpecieName)
            # else:
            #     hitTitleToSearch = hitTitle
            ### 
                
# =============================================================================
#                 Block searching for provided key words in a hit title
# =============================================================================
            
            for searchedEntity in queryList:
                for keyword in searchedEntity.searchKeywords:        ### TODO filter overlapping and similar annotations; leave the different ones 
                    # print(f'Keyword is {keyword}')
                    if keyword in hitTitleToSearch: 
                        levelPresent = False
                        
                        foundHit = BlastHit(
                            hitTitle,
                            keyword,
                            entry.descriptions[numberOfHit].e,
                            hit.accession,
                            entry.descriptions[numberOfHit].items[0].taxid,
                            titleSpecieName,
                            hit.hsps[0].align_length,
                            hit.hsps[0].identities,
                            searchedEntity)
                        
                        if entry.query in dictOfSatisfyingProteins:
                            # print("Hi, I'm here!")
                            ######
                            # Searching for the annotation of the same subclass in a sequence
                            # XXX Currently I don't check for "physycally" overlapping assignments
                            #####
                            # hitSubclass = searchedEntity.subclass
                            
                            hitClassificationValue = getattr(searchedEntity, level)
                            
                            # print(f"Hit subclass is {hitSubclass}")
                            # print(f"Already have {dictOfSatisfyingProteins[entry.query].hits[0].classificationInfo.subclass}")
                            for i in dictOfSatisfyingProteins[entry.query].hits: # TODO just write all the subclasses that you have 
                                # if i.classificationInfo.subclass == hitSubclass:
                                if getattr(i.classificationInfo, level) == hitClassificationValue:
                                    levelPresent = True
                                    break
                            
                            #####                            
                            
                            if levelPresent:
                                continue
                            else:
                                dictOfSatisfyingProteins[entry.query].hits.append(foundHit)
                                searchedEntity.incrementCounter()
                                

                        else:
                            dictOfSatisfyingProteins[entry.query] = proteinRecord(
                                entry.query,
                                foundHit
                                )
                            searchedEntity.incrementCounter()
                            
                        matchFound = True
                        break               ### XXX Now it breaks after the first "positive" consider stop after e-values comparison
            numberOfHit += 1
        # print(len(listOfSatisfyingProteins))
        
        if matchFound:
            # print("match found so +1")
            matchesFound += 1
        
        if counterBlastEntries % reportFrequency == 0:
            print(f"{counterBlastEntries} blast entries are scanned")
            print(f"{matchesFound} matches found to the moment")
            
        
        
    print(f"Total {counterBlastEntries} blast entries are processed\n{matchesFound} matches found")
    
    return dictOfSatisfyingProteins
    
    # for element in dictOfSatisfyingProteins.values():
    #     print(element)


def searchIPStsv(tsv, queryList, level = "subclass", reportFrequency = 20000):
    '''
    Searches for keywords in 6th and 13th positions of IPS produced tsv file

    Parameters
    ----------
    tsv : TYPE
        DESCRIPTION.
    queryList : TYPE
        DESCRIPTION.
    reportFrequency : TYPE, optional
        DESCRIPTION. The default is 20000.

    Returns
    -------
    None.

    '''
    
    relevantLevels = {"domain", "subclass", "class"}
    
    if level not in relevantLevels:
        raise ValueError ('Level must be one of the following ("domain", "subclass", "class")')
    elif level == "domain":
        level = "name"
    elif level == "class":
        level = "className"
        
    print(f"Level is {level}")
    
    dictOfSatisfyingProteins = {}
    counterIPSentries = 0
    matchesFound = 0
    lastProteinName = ""
    foundClassifications = set()
    tsvReader = csv.reader(tsv, delimiter="\t")
    
    
    for IPStsvRow in tsvReader:
        ########
        ### The fact that all records for one protein ID go together helps a lot
        ########
        
        matchFound = False
        
        proteinName = IPStsvRow[0]
        
        if proteinName != lastProteinName and lastProteinName != "":
            foundClassifications = set()
    
        try:
            hitTitle = IPStsvRow[12]
            strTosearchIn = ' '.join([IPStsvRow[i].lower() for i in [5,12]])          
        except IndexError:            
            hitTitle = IPStsvRow[5]
            strTosearchIn = hitTitle.lower()            
            
        # print(f'Searching in {strTosearchIn}')

        for searchedEntity in queryList:
            for keyword in searchedEntity.searchKeywords:        ### TODO filter overlapping and similar annotations; leave the different ones
                # print(f'Keyword is {keyword}')
                if keyword in strTosearchIn: 
                    
                    foundHit = IPShit(
                        hitTitle,
                        keyword,
                        IPStsvRow[8],
                        searchedEntity)
                    
                    # hitSubclass = searchedEntity.subclass
                    hitClassification = getattr(searchedEntity, level)
                    
                    ##########
                    ### Checking if we have such a subclass for a current protein
                    ##########
                    
                    
                    if not foundClassifications:
                        matchFound = True
                        foundClassifications.add(hitClassification)
                        dictOfSatisfyingProteins[proteinName] = proteinRecord(
                            proteinName,
                            foundHit
                            )
                        searchedEntity.incrementCounter()
                        break
                    elif hitClassification in foundClassifications:
                        break
                    else:
                        foundClassifications.add(hitClassification)
                        
                    
                    dictOfSatisfyingProteins[proteinName].hits.append(foundHit)
                    searchedEntity.incrementCounter()
                            
                                          
                    
                    break               ### XXX Break corresponds to search within the same searchedEntity.searchKeywords. Now it breaks after the first "positive" consider stop after e-values comparison
    # print(len(listOfSatisfyingProteins))
        if matchFound:
            # print("match found so +1")
            matchesFound += 1
            matchFound = False
        
        counterIPSentries += 1
        lastProteinName = proteinName
        
        if counterIPSentries % reportFrequency == 0:
            print(f"{counterIPSentries} IPS entries are scanned")
            print(f"{matchesFound} matches found to the moment")
        
    
    
    print(f"Total {counterIPSentries} IPS entries are processed\n{matchesFound} matches found")
    
    print(f"Dict length is {len(dictOfSatisfyingProteins)}")
    
    # for hitInfo in dictOfSatisfyingProteins.values():
    #     print(hitInfo)

    return dictOfSatisfyingProteins
    
    



def writeOutputFiles(dictToWrite, outputFolder, level="subclass", mode="BlastXML"):
    '''
    Iterates over all found proteins, checks their functional class and writes information into files.
    Each functional class has its own file.
    Creates corresponding file when first meets a new class.
    After dict of proteins is exhausted uses one more function to close all the open files.

    Parameters
    ----------
    dictToWrite : dictionary
        DESCRIPTION.
    outputFolder : string
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
# =============================================================================
#     Level check part is not necessary as I ensure that it has a correct value using argparse during data input.
# =============================================================================
    
    relevantLevels = {"domain", "subclass", "class"}
    
    print(f"Level is {level}")
    
    if level not in relevantLevels:
        raise ValueError (f'Level must be one of the following ("domain", "subclass", "class") but it is "{level}"')
    elif level == "domain":
        level = "name"
    elif level == "class":
        level = "className"
# =============================================================================

    
    dictOfFiles = {}
    dictOdDataframes = {}
    
    def createOutputFile(folder, levelName):
        if mode == "BlastXML":
            openedFile = open(folder+'/'+levelName.replace(" ", "_") + '_Blast' + '.tsv','w')
        else:
            openedFile = open(folder+'/'+levelName.replace(" ", "_") + '_IPS' + '.tsv','w')
        
        if mode == "BlastXML":
            headerLine = '\t'.join(['Sequence name','Keyword','e-value','identity','subclass','class','hit title','hit accession','hit TaxID'])+'\n'
        elif mode == "IPStsv":
            headerLine = '\t'.join(['Sequence name','Keyword','e-value','subclass','class','hit title'])+'\n'
        else:
            Exception(f'Mode is "{mode}", though it has to be one of "BlastXML" or "IPStsv"')
            
        openedFile.write(headerLine)
        return {levelName: openedFile}
    
    def createAndWriteShortTable(dataframe, levelName):
        
        def strip_seqID(string):
            match = re.search("\|.*\|", string)
            if match:
                return string[match.end()+1:]
            else:
                return string
            
        dataframe["hit title"] = dataframe["hit title"].apply(strip_seqID)
        
        if mode == "BlastXML":
            pathToWrite = outputFolder + "/" + levelName.replace(" ", "_") + "_Blast_short" + ".tsv"
            dataframe[["identity"]] = dataframe[["identity"]].astype('float')
            dataframe[["hit title", "Keyword", "e-value", "identity"]].sort_values(by=["hit title", "identity"], ascending=[True, False]).drop_duplicates(subset=["hit title"]).to_csv(pathToWrite, sep="\t", index=False)
        else:
            pathToWrite = outputFolder + "/" + levelName.replace(" ", "_") + "_IPS_short" + ".tsv"
            dataframe[["e-value"]] = dataframe[["e-value"]].astype('float')
            dataframe[["hit title", "Keyword", "e-value"]].sort_values(by=["hit title", "e-value"], ascending=[True, True]).drop_duplicates(subset=["hit title"]).to_csv(pathToWrite, sep="\t", index=False)
            
    def convertDashes(string):
        if string == "-":
            return "1e-3"
        else:
            return string
        
    
    def closeOutputFiles(filesDict):
        for file in filesDict.values():
            file.close()
    
    for protein in dictToWrite.values():
        for hit in protein.hits:
            # hitSubclass = hit.classificationInfo.subclass
            hitClassificationLevel = getattr(hit.classificationInfo, level)
            
            
            if mode == "BlastXML":
                blastString = '\t'.join([protein.name,
                                    hit.keyword,
                                    str(hit.eValue),
                                    str(hit.identity),
                                    hit.classificationInfo.subclass,
                                    hit.classificationInfo.className,
                                    hit.hitName,
                                    hit.accession,
                                    str(hit.hitTaxNum)])+'\n'
            else:
                IPSstring = '\t'.join([protein.name,
                                    hit.keyword,
                                    str(hit.eValue),
                                    hit.classificationInfo.subclass,
                                    hit.classificationInfo.className,
                                    hit.hitName])+'\n'
            
            if hitClassificationLevel in dictOfFiles:
                if mode == "BlastXML":
                    dictOfFiles[hitClassificationLevel].write(blastString)
                    dictOdDataframes[hitClassificationLevel] = dictOdDataframes[hitClassificationLevel].append(pd.DataFrame([[hit.hitName,
                                                                    hit.keyword,
                                                                    hit.eValue,
                                                                    hit.identity]],
                                                                    columns=["hit title", "Keyword", "e-value", "identity"]),
                                                                    ignore_index=True)
                    
                elif mode == "IPStsv":
                    dictOfFiles[hitClassificationLevel].write(IPSstring)
                    dictOdDataframes[hitClassificationLevel] = dictOdDataframes[hitClassificationLevel].append(pd.DataFrame([[hit.hitName,
                                                                    hit.keyword,
                                                                    convertDashes(hit.eValue)]],
                                                                    columns=["hit title", "Keyword", "e-value"]),
                                                                    ignore_index=True)
                    
            else:
                dictOfFiles.update(createOutputFile(outputFolder, hitClassificationLevel))
                if mode == "BlastXML":
                    dictOfFiles[hitClassificationLevel].write(blastString)
                    dictOdDataframes.update({hitClassificationLevel:pd.DataFrame({"hit title": [hit.hitName],
                                                                                  "Keyword": [hit.keyword],
                                                                                  "e-value": [hit.eValue],
                                                                                  "identity": [hit.identity]})})
                    
                elif mode == "IPStsv":
                    dictOfFiles[hitClassificationLevel].write(IPSstring)
                    dictOdDataframes.update({hitClassificationLevel:pd.DataFrame({"hit title": [hit.hitName],
                                                                                  "Keyword": [hit.keyword],
                                                                                  "e-value": [convertDashes(hit.eValue)]})})
    closeOutputFiles(dictOfFiles)

    for levelName, data in dictOdDataframes.items():
        createAndWriteShortTable(data, levelName)

# =============================================================================
# Body
# =============================================================================


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Searches for the entries in XML file produced by Blast, needs tsv file with specifiations what to search."
    )

    parser.add_argument(
        "-i",
        "--input_file",
        help="Blast-produced XML 2.0",
        type=argparse.FileType("r"),
        required=True,
    )
    parser.add_argument(
        "--listOfProteins",
        help="tsv file which enlists what to search in XML. Structure is as follows: \
                        [Full name	Primary name	Alias	Subclass	Class], header is present",
        type=argparse.FileType("r"),
        required=True,
    )
        
    parser.add_argument(
        "-o",
        "--outputFolder",
        help="Path to a folder where the script will create output files",
        type=str,
        default='./'
        )
    parser.add_argument(
        "-m",
        "--mode",
        help="Mode which can be 'IPStsv' or BlastXML",
        type=str,
        choices={"IPStsv", "BlastXML"},
        default="BlastXML")
    parser.add_argument(
        "-l",
        "--level",
        help='Annotation level to aggregate found molecules. Can be one of ("domain", "subclass", "class")',
        type=str,
        choices={"domain", "subclass", "class"},
        default="subclass")

    args = parser.parse_args()

    inputFile = args.input_file
    proteinTSV = args.listOfProteins
    outputFolder = args.outputFolder
    mode = args.mode
    aggregationLevel = args.level

    proteinsToSearch = readProteinTSV(proteinTSV)
    
    if mode == "BlastXML":

        dictOfFoundProteins = searchBlastXML(inputFile,
                                             proteinsToSearch,
                                             level = aggregationLevel)
        
        writeOutputFiles(dictOfFoundProteins, outputFolder, level = aggregationLevel)
    
    elif mode == "IPStsv":
        dictOfFoundProteins = searchIPStsv(inputFile,
                                           proteinsToSearch, 
                                           level = aggregationLevel)
        writeOutputFiles(dictOfFoundProteins, outputFolder,level = aggregationLevel, mode = "IPStsv")
    
    inputFile.close()

    var1 = "kappaPride"

    print("It's okay")

    # dir(__builtins__)

    # print(globals())

    # print(locals())


