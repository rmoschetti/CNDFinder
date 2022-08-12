__version__ = '1.9'
__all__ = [
    'CndFinder', 'CPhiFinder'
]

import logging, sys
from sage.all import *

Log_Format = "%(asctime)s - %(message)s"

logging.basicConfig(stream = sys.stdout, 
                    format = Log_Format,
                    datefmt = '%I:%M:%S',
                    level = logging.DEBUG)

logger = logging.getLogger()
logging.getLogger().setLevel(logging.DEBUG)

############################################################################
#The main function which computes the combinatorial non-degeneracy invariant
############################################################################

def CndFinder(IntersectionMatrix,GeneratorsList, *args, **kwargs):
    #CndFinder(IntersectionMatrix,GeneratorsList) computes the combinatorial non-degeneracy invariant cnd(S,R) of the system of curves R
    #whose intersection matrix is IntersectionMatrix, given a list of generators of Num(S) in GeneratorsList. 
    #input: IntersectionMatrix: the intersection matrix of the collection of smooth rational curves R
    #       GeneratorsList: a list of arrays, each array is the coordinate vector of an element of a basis of Num(S), 
    #                       written as a combination of curves in R
    #       SubGraphsToSearch:={'A':[...], 'D':[...], 'E':[...]} optional 
    #                       (default: {'A':[1,2,3,4,5,6,7,8],'D':[4,5,6,7,8],'E':[6,7,8]}). 
    #                       A dictionary specifying the maximum rank of extended Dynkin diagrams ADE the program should look for
    #       ShiftPrint: optional (default: 1). Shift the labelling of the curves by this integer in the report. 
    #                      This does not affect the value for cnd(S,R).
    #output: a dictionary containing the following data. 
    #       'EllipticFibrations' contains the list of elliptic configuration grouped in elliptic fibrations.
    #       'Cnd', containing the value cnd(S,R). 
    #       'SaturatedSequences', a list containing all the R-saturated sequences of vectors in HF(S,R) found by the program.
    
    logger.info("*---- CndFinder starts. ----*")
       
    #Get the optional arguments
    SubGraphsToSearch = kwargs.get('SubGraphsToSearch', {'A':[1,2,3,4,5,6,7,8],'D':[4,5,6,7,8],'E':[6,7,8]})
    ShiftPrint = kwargs.get('ShiftPrint', 1)
        
    ###############################################################################################
    ## Begin of STEP 1
    ###############################################################################################
    
    logger.info("*---- Begin of STEP 1 - Find all the elliptic configurations supported on the given curves ---- *")
    
    #Generate the intersection matrices of the diagrams specified in SubGraphsToSearch. 
    #The diagrams are dual graphs to elliptic configurations.
    
    TemplateDiagrams=[]
    for item in SubGraphsToSearch['A']:
        TemplateDiagrams.append(extendedDynkinGraph('A',item))
    for item in SubGraphsToSearch['D']:
        TemplateDiagrams.append(extendedDynkinGraph('D',item))
    for item in SubGraphsToSearch['E']:
        TemplateDiagrams.append(extendedDynkinGraph('E',item))
        

    #TempEllipticConfiguration[i] will be in turn a list of all subgraphs of type 'TemplateDiagrams[i]', which support an elliptic configuration.
    TempEllipticConfigurations=[] 

    for item in TemplateDiagrams:        
        #Look for all the elliptic configurations whose dual graph is item['Matrix'], for item in TemplateDiagrams
        logger.info("\tFind elliptic configurations of type " + item['Type'])
        TempEllipticConfigurations.append(FindSubgraphs(IntersectionMatrix,item['Matrix'],[[]]))

        
    logger.info("*---- End of STEP 1. ----*")    
    
    ###############################################################################################
    ## End of STEP 1
    ###############################################################################################
    
    ###############################################################################################
    ## Begin of STEP 2
    ###############################################################################################
    
    logger.info("*---- Begin of STEP 2 - Check fibers and half-fibers ----*")
    
    #Initialize the list which will contain the elliptic configurations
    EllipticConfigurationList=[]
    
    for index,item in enumerate(TempEllipticConfigurations):
        #Subdivide in fibers and half-fibers the elliptic configurations in item:=TempEllipticConfigurations[index]
        #we use the fact that the list TemplateDiagrams is indexed in the same way as TempEllipticConfigurations, so 
        #item:=TempEllipticConfigurations[index] contains the list of elliptic configuration of type TemplateDiagrams[index]
        tempFiberHalfFiber=catFibersHalfFibers(item,TemplateDiagrams[index],IntersectionMatrix,GeneratorsList)
        #Given item:=TempEllipticConfigurations[index], tempFiberHalfFiber will be a list [data of the fibers, data of the half-fibers]
        #Check the function TempEllipticConfigurations for more information.

        for item in tempFiberHalfFiber:
            logger.info("\tType: " + item['Type'] + ". " + 
                        "Number: " + str(len(item['list']))+ ". ")
            
        #Add the elliptic configurations to the main list. 
        EllipticConfigurationList.extend(tempFiberHalfFiber)

    logger.info("*---- End of STEP 2. ----*")    
    
    ###############################################################################################
    ## End of STEP 2
    ###############################################################################################
    
    ###############################################################################################
    ## Begin of STEP 3
    ###############################################################################################

    logger.info("*---- Begin of STEP 3 - Collect elliptic configurations in elliptic fibrations. ----*")

    #Initialize the dictionary of elliptic fibrations, 
    #which will group together equivalent elliptic configurations of EllipticConfigurationList
    EllipticFibrationsList={}

    #Cycle through all the elliptic configurations of EllipticConfigurationList, 
    #and group together the ones which intersect giving zero: hence belong to the same elliptic fibration.
    #The list EllipticConfigurationList contains all the elliptic configurations divided by their type. 
    #So first cycle through each possible type of elliptic configuration, named by "item".
    for item in EllipticConfigurationList:
        
        #Then, consider each elliptic configuration of the specific type, which is listed under item['list']
        for ec in item['list']:
            
            #Check if ec is already contained in an elliptic fibration found previously 
            IsContained=False
            #The instruction sum([...],[]) just flattens the nested list of representatives of elliptic fibrations
            for v in sum([EllipticFibrationsList[s]['list'] for s in EllipticFibrationsList],[]):
                #ec is already in the list if and only if there already is a representative v which intersects ec in 0
                if intProduct(v,ec,IntersectionMatrix)==0:
                    IsContained=True
                    break
            
            
            if (IsContained==False):
                #In case ec is new, find all the elliptic configurations disjoint from ec
                DisjointEllipticConfigurationFound=GetDisjointEllipticConfiguration(ec,EllipticConfigurationList,IntersectionMatrix)


                #For the final report, store in TempFibersVector the list all the elliptic configurations in this elliptic fibration
                #expressed in vector form
                TempFibersVector=[GetEllipticConfiguration(EllipticConfigurationList,s) for s in DisjointEllipticConfigurationFound]

                #Finalize the data of the elliptic fibration to be stored. 
                #Note that only 'Representative' will be used again in the computation of Cnd, 
                #the elliptic configurations equivalent to 'Representative' are stored in 'FibersVector' and 'FibersSupport'.
                FibrationDataToAdd={'Representative':ec,
                                'FibersVector':TempFibersVector,
                                'FibersSupport':ListVectorToListSupport([TempFibersVector],ShiftPrint)
                               }

                #Generate a name to index the elliptic fibration, according to the type of its singular fibers
                NameFibration="(" + sequenceToDescriptionString(DisjointEllipticConfigurationFound,EllipticConfigurationList, " ", " + ") + ")"
            
                #Add the data of the fibration to the list "EllipticFibrationsList"
                if (NameFibration in EllipticFibrationsList):
                    EllipticFibrationsList[NameFibration]['FibrationsData'].append(FibrationDataToAdd)
                    EllipticFibrationsList[NameFibration]['list'].append(ec)
                else:
                    logger.info("\t Found new elliptic configuration type: " + NameFibration)
                    EllipticFibrationsList[NameFibration]={'Type': NameFibration, 
                                                            'RepresentativesType': item['Type'],
                                                             #we put under the name 'list', the list of representatives
                                                            'list': [ec],
                                                            'FibrationsData':[FibrationDataToAdd]}
          

    #Convert the dictionary "EllipticFibrationsList" in a list
    EllipticFibrationsList=list(EllipticFibrationsList.values())
    
    logger.info("*---- End of STEP 3. ----*")    
    
    ###############################################################################################
    ## End of STEP 3
    ###############################################################################################
    
    ###############################################################################################
    ## Begin of STEP 4
    ###############################################################################################
    
    logger.info("*---- Begin of STEP 4 - Compute the longest non-degenerate isotropic sequence with fibrations of only one type ----*")
    
    #Compute the longest non-degenerate isotropic sequence containing only elements of a given fibration type, i.e. each element of EllipticFibrationsList. This step, while not necessary, helps speeding up the computations.
    
    for i in range(0,len(EllipticFibrationsList)):
        TempListOfSequences=[[]]     
        
        #Try to add an elliptic fibration of type i to the sequences in TempListOfSequences, for up to 10 times.
        for maxNumber in range(0,11):
            TempListOfSequences=AddRepresentativeListToSequenceList(TempListOfSequences,i,EllipticFibrationsList,IntersectionMatrix)

            #Stop when the list of resulting sequences has length 0 (so it is not possible to add more fibrations of Type i)
            if (len(TempListOfSequences)==0):
                break 
        
        #Save this as the maximum number of elliptic fibrations having Type i which can appear together in an isotropic sequence.
        EllipticFibrationsList[i]['max']=maxNumber
                   
        logger.info("\tType: " + EllipticFibrationsList[i]['Type'] + 
                    "\n\t\tNumber: " + str(len(EllipticFibrationsList[i]['list'])) + "; " +
                    "\n\t\tMax: " + str(EllipticFibrationsList[i]['max']) + "; ")
    
    logger.info("*---- End of STEP 4. ----*")    
    
    ###############################################################################################
    ## End of STEP 4
    ###############################################################################################
    
    ###############################################################################################
    ## Begin of STEP 5
    ###############################################################################################
    
    logger.info("*---- Begin of STEP 5 - Compute the cnd-invariant ----*")
    
    #Initialize the list which will collect the final results
    ResultList=[]
    
    #Call the recursive function which will compute the cnd-invariant, starting to extend the empty sequence [[]] with the zero-th element in EllipticFibrationsList and return the final results
    InitialEllipticFibrationType=0

    CndFound=RecursiveCnd([[]],InitialEllipticFibrationType,ResultList,EllipticFibrationsList,IntersectionMatrix,ShiftPrint)
    
    logger.info("*---- End of STEP 5. ----*")    
    
    ###############################################################################################
    ## End of STEP 5
    ###############################################################################################
    
    ###############################################################################################
    ## Begin of STEP 6
    ###############################################################################################
    
    logger.info("*---- Begin of STEP 6 - Find the R-saturated sequences ----*")
    
    #Initialize the list which will contain the R-saturated sequences
    SaturatedSequences=[]
    

    #Scan all the results, and save only the R-saturated isotropic sequences
    for item in ResultList:

        #Create a list of candidate isotropic sequences which may extend isotropic sequences in item. 
        #The instruction sum([...],[]) just flattens the list resulting from the calculation [...]
        PossibleExtensions=sum([s['ListSequencesSupport'] for s in ResultList if IsStrictlyContained(item['TypeVector'],s['TypeVector'])],[])

        #Initialize some temporary list where to store the R-saturated sequences in item
        ListSaturatedSequencesSupportNew=[]
        ListSaturatedSequencesVectorNew=[]

        #Scan all the sequences of type "Type"
        for i in range(0,len(item['ListSequencesSupport'])):
            #Check if they are R-saturated, i.e they are not subsequences of any elements in PossibleExtensions
            if not (IsInList(item['ListSequencesSupport'][i],PossibleExtensions)==True):
                ListSaturatedSequencesSupportNew.append(item['ListSequencesSupport'][i])
                ListSaturatedSequencesVectorNew.append(item['ListSequencesVector'][i])

        if len(ListSaturatedSequencesSupportNew)>0:
            SaturatedSequences.append({'Found cnd':item['Found cnd'],
                           'Type':item['Type'],
                           'Cardinality': len(ListSaturatedSequencesVectorNew),
                           'ListSequencesVector': ListSaturatedSequencesVectorNew,
                           'ListSequencesSupport': ListSaturatedSequencesSupportNew,
                           'TypeVector': item['TypeVector']})
    
    
    
    logger.info("*---- End of STEP 6. ----*")    
    
    ###############################################################################################
    ## End of STEP 6
    ###############################################################################################
 
    logger.info("*---- Summary of the results. ----*")
    logger.info("*---- This is the combinatorial cnd: " + str(CndFound) + " ----*")
    logger.info("*---- Here is a list of the saturated sequences found ----*")
    for item in SaturatedSequences:
        logger.debug("\tLength: " + str(item['Found cnd']))
        logger.debug("\tType: " + str(item['Type']))
        logger.debug("\tCardinality: "  + str(item['Cardinality']))
        logger.debug("\tExample: "  + str(item['ListSequencesSupport'][0]))
        logger.debug("------------------------------------------------------------------------")
    
    #Initialize the dictionary which will contain the final data, wrap them up and return them to the user
    FinalData={}
    
    FinalData['Cnd']=CndFound
    FinalData['SaturatedSequences']=SaturatedSequences
    FinalData['EllipticFibrations']=EllipticFibrationsList
    
    NicePrint(FinalData,GeneratorsList,IntersectionMatrix)

    return FinalData
    

    
############################################################################
#The main function that computes the combinatorial Phi-invariant
############################################################################

def CPhiFinder(IntersectionMatrix,GeneratorsList, DivisorList, *args, **kwargs):
    #CPhiFinder(IntersectionMatrix,GeneratorsList) computes the combinatorial Phi-invariant Phi(S,R,D) of a divisor D on S with
    #respect to the smooth rational curves in R, whose intersection matrix is IntersectionMatrix, 
    #given a list of generators of Num(S) in GeneratorsList. 
    #input: IntersectionMatrix: the intersection matrix of the collection of smooth rational curves in R
    #       GeneratorsList: a list of arrays, each array is the coordinate vector of an element of a basis of Num(S), 
    #                       written as a combination of curves in R
    #       SubGraphsToSearch:={'A':[...], 'D':[...], 'E':[...]} optional 
    #                       (default: {'A':[1,2,3,4,5,6,7,8],'D':[4,5,6,7,8],'E':[6,7,8]}). 
    #                       A dictionary specifying the maximum rank of extended Dynkin diagrams ADE the program should look for
    #       DivisorList: a list of divisors of which we want to compute the combinatorial Phi-invariants.
    #                           Each divisor is represented by a list of rational numbers which are the coefficients of a linear combination of the 
    #                           smooth rational curves in R. Note that the program does not check if such linear combination is in Num(S) or if it is big.    
    #                              
              
    
    logger.info("*---- CPhiFinder starts. ----*")
       
    #Get the optional arguments
    SubGraphsToSearch = kwargs.get('SubGraphsToSearch', {'A':[1,2,3,4,5,6,7,8],'D':[4,5,6,7,8],'E':[6,7,8]})
    ShiftPrint = kwargs.get('ShiftPrint', 1)
        
    ###############################################################################################
    ## Begin of STEP 1 - identical to STEP 1 of CndFinder
    ###############################################################################################
    
    logger.info("*---- Begin of STEP 1 - Find all the elliptic configurations supported on the given curves ---- *")
    
    #Generate the intersection matrices of the diagrams specified in SubGraphsToSearch. 
    #The diagrams are dual graphs to elliptic configurations.
    
    TemplateDiagrams=[]
    for item in SubGraphsToSearch['A']:
        TemplateDiagrams.append(extendedDynkinGraph('A',item))
    for item in SubGraphsToSearch['D']:
        TemplateDiagrams.append(extendedDynkinGraph('D',item))
    for item in SubGraphsToSearch['E']:
        TemplateDiagrams.append(extendedDynkinGraph('E',item))
        

    #TempEllipticConfiguration[i] will be in turn a list of all subgraphs of type 'TemplateDiagrams[i]', which support an elliptic configuration.
    TempEllipticConfigurations=[] 

    for item in TemplateDiagrams:        
        #Look for all the elliptic configurations whose dual graph is item['Matrix'], for item in TemplateDiagrams
        logger.info("\tFind elliptic configurations of type " + item['Type'])
        TempEllipticConfigurations.append(FindSubgraphs(IntersectionMatrix,item['Matrix'],[[]]))

        
    logger.info("*---- End of STEP 1. ----*")    
    
    ###############################################################################################
    ## End of STEP 1
    ###############################################################################################
    
    ###############################################################################################
    ## Begin of STEP 2 - identical to STEP 2 of CndFinder
    ###############################################################################################
    
    logger.info("*---- Begin of STEP 2 - Check the half-fibers ----*")
    
    #Initialize the list which will contain the elliptic configurations
    EllipticConfigurationList=[]
    
    for index,item in enumerate(TempEllipticConfigurations):
        #Subdivide in fibers and half-fibers the elliptic configurations in item:=TempEllipticConfigurations[index]
        #we use the fact that the list TemplateDiagrams is indexed in the same way as TempEllipticConfigurations, so 
        #item:=TempEllipticConfigurations[index] contains the list of elliptic configuration of type TemplateDiagrams[index]
        tempFiberHalfFiber=catFibersHalfFibers(item,TemplateDiagrams[index],IntersectionMatrix,GeneratorsList)
        #Given item:=TempEllipticConfigurations[index], tempFiberHalfFiber will be a list [data of the fibers, data of the half-fibers]
        #Check the function TempEllipticConfigurations for more information.

        for item in tempFiberHalfFiber:
            logger.info("\tType: " + item['Type'] + ". " + 
                        "Number: " + str(len(item['list']))+ ". ")
            
        #The function catFibersHalfFibers returns both the fibers and the half-fibers.
        #The coefficients of the fibers have already been divided by two in catFibersHalfFibers, hence, for each fiber F we already have
        #the coefficients of a half-fiber H such that 2H=F
        EllipticConfigurationList.extend(tempFiberHalfFiber)

    logger.info("*---- End of STEP 2. ----*")    
    
    ###############################################################################################
    ## End of STEP 2
    ###############################################################################################
    
    ###############################################################################################
    ## Begin of STEP 3 - identical to STEP 3 of CndFinder
    ###############################################################################################

    logger.info("*---- Begin of STEP 3 - Collect elliptic configurations in elliptic fibrations. ----*")

    #Initialize the dictionary of elliptic fibrations, 
    #which will group together equivalent elliptic configurations of EllipticConfigurationList
    EllipticFibrationsList={}

    #Cycle through all the elliptic configurations of EllipticConfigurationList, 
    #and group together the ones which intersect giving zero: hence belong to the same elliptic fibration.
    #The list EllipticConfigurationList contains all the elliptic configurations divided by their type. 
    #So first cycle through each possible type of elliptic configuration, named by "item".
    for item in EllipticConfigurationList:
        
        #Then, consider each elliptic configuration of the specific type, which is listed under item['list']
        for ec in item['list']:
            
            #Check if ec is already contained in an elliptic fibration found previously 
            IsContained=False
            #The instruction sum([...],[]) just flattens the nested list of representatives of elliptic fibrations
            for v in sum([EllipticFibrationsList[s]['list'] for s in EllipticFibrationsList],[]):
                #ec is already in the list if and only if there already is a representative v which intersects ec in 0
                if intProduct(v,ec,IntersectionMatrix)==0:
                    IsContained=True
                    break
            
            
            if (IsContained==False):
                #In case ec is new, find all the elliptic configurations disjoint from ec
                DisjointEllipticConfigurationFound=GetDisjointEllipticConfiguration(ec,EllipticConfigurationList,IntersectionMatrix)


                #For the final report, store in TempFibersVector the list all the elliptic configurations in this elliptic fibration
                #expressed in vector form
                TempFibersVector=[GetEllipticConfiguration(EllipticConfigurationList,s) for s in DisjointEllipticConfigurationFound]

                #Finalize the data of the elliptic fibration to be stored. 
                #Note that only 'Representative' will be used again in the computation of Cnd, 
                #the elliptic configurations equivalent to 'Representative' are stored in 'FibersVector' and 'FibersSupport'.
                FibrationDataToAdd={'Representative':ec,
                                'FibersVector':TempFibersVector,
                                'FibersSupport':ListVectorToListSupport([TempFibersVector],ShiftPrint)
                               }

                #Generate a name to index the elliptic fibration, according to the type of its singular fibers
                NameFibration="(" + sequenceToDescriptionString(DisjointEllipticConfigurationFound,EllipticConfigurationList, " ", " + ") + ")"
            
                #Add the data of the fibration to the list "EllipticFibrationsList"
                if (NameFibration in EllipticFibrationsList):
                    EllipticFibrationsList[NameFibration]['FibrationsData'].append(FibrationDataToAdd)
                    EllipticFibrationsList[NameFibration]['list'].append(ec)
                else:
                    logger.info("\t Found new elliptic configuration type: " + NameFibration)
                    EllipticFibrationsList[NameFibration]={'Type': NameFibration, 
                                                            'RepresentativesType': item['Type'],
                                                             #we put under the name 'list', the list of representatives
                                                            'list': [ec],
                                                            'FibrationsData':[FibrationDataToAdd]}
          

    #Convert the dictionary "EllipticFibrationsList" in a list
    EllipticFibrationsList=list(EllipticFibrationsList.values())
    
    logger.info("*---- End of STEP 3. ----*")    
    
    ###############################################################################################
    ## End of STEP 3
    ###############################################################################################
    
    ###############################################################################################
    ## Begin of FINAL STEP
    ###############################################################################################
    
    logger.info("*---- Begin of FINAL STEP - Find the combinatorial Phi-invariant ----*")
    
    #In the list Result we store the combinatorial Phi-invariants of the divisors in DivisorList
    Result=[]
    
    for Divisor in DivisorList:
        logger.debug("\tComputing the combinatorial Phi-invariant of " + str(Divisor))
        #The starting value for the variable PhiInvariant is set to be infinity and it is replaced by a finite value after the first iteration of the cycle below
        PhiInvariant=infinity

        for item in EllipticFibrationsList:
            for t in item['list']:
                #Compute the product of Divisor with the Half-Fiber t
                tempProduct=intProduct(t,Divisor,IntersectionMatrix)

                if (tempProduct<=PhiInvariant):
                    PhiInvariant=tempProduct
                    logger.debug("\t\tNew lower bound found: "  + str(PhiInvariant) + 
                                 ". Given by the half-fiber of type " + item['Type'] + 
                                 " with coefficients " + str(t))

        logger.info("\t\tValue of the Phi-invariant found: "+ str(PhiInvariant))
        Result.append({'Divisor':Divisor,'CPhi':PhiInvariant})
    logger.info("*---- End of FINAL STEP. ----*")
    
    return Result
    
    ###############################################################################################
    ## End of FINAL STEP
    ###############################################################################################
    
    
##########################################################
#Auxiliary functions called in various steps
##########################################################

    
def intProduct(P,Q,M):
    #intProduct(P,Q,M) computes the internal product defined by the matrix M of the vectors given by P and Q.
    #Input: square matrix M of dimension n, lists P, Q of length n.
    #Output: a number
    return matrix(P)*M*transpose(matrix(Q))  


def GetEllipticConfiguration(Diagrams,d):
    #GetEllipticConfiguration takes a pair of indices d, which identify an element in the list Diagrams, and return the same element explicitly written as a vector
    
    return Diagrams[d[0]]['list'][d[1]]


    

##########################################################
# Auxiliary functions called in STEP 1
##########################################################

def extendedDynkinGraph(Type,Order):
    #ExtendedDynkinGraph generates the intersection matrix of an extended Dynkin diagram of Type "A", "D" or "E" and of specific Order.
    
    Result={}
    Result['Type']=Type + str(Order)
    
    if Type=='A':
        Result['Weights']=[1]*(Order+1)
        
        if (Order==1):
            Result['Matrix']=matrix([[-2,2],[2,-2]])
            
        else:
            Result['Matrix']=(-2)*identity_matrix(Order+1)
            Result['Matrix'][0,-1]=1
            Result['Matrix'][-1,0]=1
            for i in range(0,Order):
                Result['Matrix'][i+1,i]=1
                Result['Matrix'][i,i+1]=1
        
    if Type=='D':
        Result['Weights']=[2]*(Order+1)
        Result['Weights'][0]=1
        Result['Weights'][1]=1
        Result['Weights'][2]=1
        Result['Weights'][3]=1
        
        Result['Matrix']=(-2)*identity_matrix(Order+1)
        Result['Matrix'][0,4]=1
        Result['Matrix'][1,4]=1
        Result['Matrix'][4,0]=1
        Result['Matrix'][4,1]=1
        Result['Matrix'][2,-1]=1
        Result['Matrix'][3,-1]=1
        Result['Matrix'][-1,2]=1
        Result['Matrix'][-1,3]=1
        for i in range(4,Order):
            Result['Matrix'][i+1,i]=1
            Result['Matrix'][i,i+1]=1
        
    if Type=='E':
        if (Order==6):
            Result['Weights']=[1,2,3,2,1,2,1]
            Result['Matrix']=matrix([
                            [-2, 1, 0, 0, 0, 0, 0],
                            [1, -2, 1, 0, 0, 0, 0],
                            [0, 1, -2, 1, 0, 1, 0],
                            [0, 0, 1, -2, 1, 0, 0],
                            [0, 0, 0, 1, -2, 0, 0],
                            [0, 0, 1, 0, 0, -2, 1],
                            [0, 0, 0, 0, 0, 1, -2]
                        ])   
        if (Order==7):
            Result['Weights']=[2,3,4,3,2,1,2,1]
            Result['Matrix']=matrix([
                            [-2, 1, 0, 0, 0, 0, 0, 1],
                            [1, -2, 1, 0, 0, 0, 0, 0],
                            [0, 1, -2, 1, 0, 0, 1, 0],
                            [0, 0, 1, -2, 1, 0, 0, 0],
                            [0, 0, 0, 1, -2, 1, 0, 0],
                            [0, 0, 0, 0, 1, -2, 0, 0],
                            [0, 0, 1, 0, 0, 0, -2, 0],
                            [1, 0, 0, 0, 0, 0, 0, -2]
                        ])
        if (Order==8):
            Result['Weights']=[2,4,6,5,4,3,2,3,1]
            Result['Matrix']=matrix([
                            [-2, 1, 0, 0, 0, 0, 0, 0, 0],
                            [1, -2, 1, 0, 0, 0, 0, 0, 0],
                            [0, 1, -2, 1, 0, 0, 0, 1, 0],
                            [0, 0, 1, -2, 1, 0, 0, 0, 0],
                            [0, 0, 0, 1, -2, 1, 0, 0, 0],
                            [0, 0, 0, 0, 1, -2, 1, 0, 0],
                            [0, 0, 0, 0, 0, 1, -2, 0, 1],
                            [0, 0, 1, 0, 0, 0, 0, -2, 0],
                            [0, 0, 0, 0, 0, 0, 1, 0, -2]
                        ])
    return Result






def FindSubgraphs(M,N,ListSubgraphsFound):
    #FindSubgraphs takes two square matrices M=(a_{i,j}) and N, and find all the ordered subset of integers {i_1, ... i_k} such that the matrix a_{i,j} with i,j \in {i_1, ... i_k} is equal to N.
    #If M is the intersection matrix of a graph G, this is equivalent to finding all the subgraphs of G with intersection matrix with respect to G equal to N
    #It does this by recursion keeping the temporary results in ListSubgraphsFound, and finally returns the results as a list. 
    
    #If at a certain point in the recursion the list ListSubgraphsFound is empty, then it means that M has no submatrix isomorphic to N in the above sense, so the function returns []
    if len(ListSubgraphsFound)==0:
        return []
    
    #If the ListSubgraphsFound contains lists with as many indexes as the size of N, the search is complete, and the results can be returned. Otherwise call recursively the function FindSubgraphs
    if (len(ListSubgraphsFound[0])==len(N[0])):
        #Remove the duplicates in ListSubgraphsFound and return them
        Result=[]
        seen = []
        for i in ListSubgraphsFound:
            if not sorted(i) in seen:
                #Notice that the sort function is used only to find duplicates, the results are kept in the original order to retain information about the order of the vertices in N.
                seen.append(sorted(i))
                Result.append(i)
                
        return Result  

    
    #If we arrive at this point of the function, it means we still have to add some indexes to the lists in ListSubgraphsFound.
    #Initialize the list which will contain the extensions of the lists in ListSubgraphsFound with the new indexes
    NewListSubgraphsFound=[]
    
    #We now set L in {1,...,size(N)} and look at the submatrix (b_{i,j})_{0<i,j<L+1} of N. L is an integer corresponding to the index of row and column of N we are considering, i.e. in terms of graphs the new vertex of N we are using to extend element in ListSubgraphsFound.
    L=len(ListSubgraphsFound[0])
    
    for item in ListSubgraphsFound:
        #We will compile the list indexes suitable to be added to the list "item".
        #There are two conditions for an index "s" to be chosen:
        #(1) it is not already in item;
        #(2) M[item[k]][s] == N[k][L], for k=0...L-1.
        #    in terms of graphs, the vertex indexed by "s" must intersect the first k=0...L-1 chosen vertices of item (i.e. item[k]) in the same way as the vertex L-th intersects the first k vertices in N.
        #    This is summarized with the instruction (M[item,s] == N[list(range(L)),L])
        for C in [s for s in range(0,len(M[0])) if (not s in item) and (M[item,s] == N[list(range(L)),L])]:
            NewListSubgraphsFound.append(item+[C])
            
    return FindSubgraphs(M,N,NewListSubgraphsFound)                  

    

##########################################################
# Auxiliary functions called in STEP 2
##########################################################

def catFibersHalfFibers(FoundList,Template,IntersectionMatrix,GeneratorsList):
    #catFibersHalfFibers start with a list of elliptic configurations (whose support has dual graph Template) and sorts them deciding if they are fibers or half-fibers
    #input: FoundList: list of configurations to sort 
    #       Template: template of the diagram (encodes the support, and the appropriate coefficients of the elliptic configuration)
    #       IntersectionMatrix: intersection matrix
    #       GeneratorsList: list of the generators of Num(S)
    #output: A list [tempFiberList,tempHalfFiberList] containing the list of Fibers and Half-fibers expressed in vector notation
    tempFiberList=[]
    tempHalfFiberList=[]
    for subGraph in FoundList:       
    
        #Convert the subgraph to vector notation, in order to compute the intersection.
        subGraphVector=[0.0]*len(IntersectionMatrix[0])
        for i in range(0,len(subGraph)):
                subGraphVector[subGraph[i]]=Template['Weights'][i]/2


        #Check if subGraphVector is a fiber or a half-fiber
        if (isFiber(IntersectionMatrix,subGraphVector,GeneratorsList)):
            tempFiberList.append(subGraphVector)
        else:
            tempHalfFiberList.append([2*s for s in subGraphVector])
        
    Result=[{'Type': Template['Type'] + 'F','list':tempFiberList},{'Type': Template['Type'] + 'HF','list':tempHalfFiberList}]   

    #Remove the empty entries in Result (e.g. if the given Template has no half-fibers)
    Result=[s for s in Result if len(s['list'])>0]    
                
    return Result


def isFiber(M,VectorToCheck,GeneratorsList):
    #isFiber(VectorToCheck,GeneratorsList) checks if the elliptic configuration defined by VectorToCheck is a Fiber, by checking that the internal product of VectorToCheck with all the generators of GeneratorsList is an integer.
    #Input: A matrix M, for the scalar product
    #       a vector VectorToCheck defining the subgraph in terms of the matrix M
    #       a list of generators of the numerical group, in the same basis in which M is defined
    #Output: True or False, depending if VectorToCheck is a fiber (i.e. its class belongs to 2Num(S)) or not
    
    for item in GeneratorsList:
        if not intProduct(VectorToCheck,item,M)[0][0] in ZZ:
            return False
    return True

    

    
##########################################################
# Auxiliary functions called in STEP 3
##########################################################  

def GetDisjointEllipticConfiguration(IndexStart,EllipticConfigurationList,Matrix):
    #GetDisjointEllipticConfiguration(IndexStart,EllipticConfigurationList,Matrix) finds all the elliptic configurations in EllipticConfigurationList which are disjoint from IndexStart. The matrix Matrix defines the intersection product
    Result=[]
    for g in range(0,len(EllipticConfigurationList)):
        for el in range(0,len(EllipticConfigurationList[g]['list'])):
            if (intProduct(IndexStart,GetEllipticConfiguration(EllipticConfigurationList,[g,el]),Matrix) == 0):
                Result.append([g,el])
    return Result


##########################################################
# Auxiliary functions called in STEP 4 and STEP 5
##########################################################

def AddRepresentativeListToSequenceList(List,Index,Diagrams,Matrix):
    #AddRepresentativeListToSequenceList(List,Index,Diagrams,Matrix) takes a list "List" of isotropic sequences, and compile a new list of isotropic sequences by appending the elliptic configuration presents in "Diagrams", starting from the position "Index"
    NewList=[]
    
    #Determine the starting position of the elements which will be added to the isotropic sequences in List
    for item in List:
        if (len(item)==0):
            Start=0
        elif (item[-1][0] < Index):
            Start=0
        else:
            Start=item[-1][1]+1
        
        #For every element from the position "Start" on, check if it is compatible (i.e., it can be appended to the isotropic sequence "item" to create a longer isotropic sequence). If it is, append it to "item" and put the new sequence in "NewList".
        for i in range(Start,len(Diagrams[Index]['list'])):
            if (isCompatible([Index,i],item,Diagrams,Matrix)):
                NewList.append(item+[[Index,i]])
    
    return NewList


def isCompatible(newd,sequence,Diagrams,Matrix):
    #isCompatible returns True if the curve indexed by "newd" in Diagrams can be added to the isotropic sequence "sequence" (which happens when the product of the curve indexed by newd with all the other curves in sequence is 1), False otherwise. 
    
    for item in sequence:
        if intProduct(GetEllipticConfiguration(Diagrams,item),GetEllipticConfiguration(Diagrams,newd),Matrix)!=1:
            return False
    return True

##########################################################
# Auxiliary functions called in STEP 5
##########################################################

def RecursiveCnd(OldList,CurrentSubgraph,ResultList,EllipticConfigurationList,IntersectionMatrix,ShiftPrint):
    #RecursiveCnd(...) computes cnd(S,R), where R is a set of curves of S which intersect as IntersectionMatrix.
    #The strategy consists in adding recursively elements to the list OldList. 
    #The arguments OldList, CurrentSubgraph are used to trace the progress of the recursion. 
    #The arguments ResultList are used for storing the results for the final report
    #The arguments EllipticConfigurationList,IntersectionMatrix,ShiftPrint are the fixed data needed for the computations

    
    #Initialize ListNewCND, with the current lower bound for cnd. Adding len(OldList[0]) is redundant.
    ListNewCND=[len(OldList[0])] #
    
    #Consider all the remaining type of elliptic fibrations, where i is the index in the array of elliptic fibrations type EllipticConfigurationList
    for i in range(CurrentSubgraph,len(EllipticConfigurationList)):   
        
        #CurrentNumber represents how many elliptic fibrations of type CurrentSubgraph are currently in the list OldList. This has to be set to 1 when i is strictly greater than CurrentSubgraph, because in this case i is a new type of elliptic fibration. This is used only to do a check against EllipticConfigurationsList[i]['max'], i.e. the maximum number of diagrams of type i in an isotropic sequence, computed in Step 4.
        CurrentNumber=len([s for s in OldList[0] if s[0]==i])+1
       
            
        logger.debug("Current sequence: " + sequenceToDescriptionString(OldList[0],EllipticConfigurationList," x ", ", ") + 
                    "Cardinality: " + str(len(OldList)) +  
                    ".\n\t\t\t\t I am trying to add: " + EllipticConfigurationList[i]['Type'] + 
                    ".\n\t\t\t\t After the addition, there will be " + str(CurrentNumber) + " "
                    + EllipticConfigurationList[i]['Type'] +" in the list.")
        
        #Check if we can skip this case, because it puts more elliptic configurations in a sequence than the maximum possible
        if (CurrentNumber>EllipticConfigurationList[i]['max']):
            logger.debug("\t\t\t Maximum number of "+ EllipticConfigurationList[i]['Type'] + " exceeded. Moving on.")
            logger.debug("------------------------------------------------------------------------")
            continue
        
        
        #Add to the current sequences in OldList all the (representative of the) elliptic fibrations indexed by i in the array EllipticConfigurationList.
        NewList=AddRepresentativeListToSequenceList(OldList,i,EllipticConfigurationList,IntersectionMatrix)
        
        #If the previous step succeeded, NewList is non-empty. In this case, there are ways to increase the length of the sequence. So we can log them and call again the recursion.
        if (len(NewList)>0):
            #Add the new length to the list of cnd found
            ListNewCND.append(len(NewList[0]))
            logger.debug("Found cnd=" + str(len(NewList[0])))
            logger.debug("\tType=" + sequenceToDescriptionString(NewList[0],EllipticConfigurationList," x ", ", "))
            logger.debug("\tCardinality=" + str(len(NewList)))
            logger.debug("\tExample=" + sequenceToVerticesListToString(NewList[0],EllipticConfigurationList,ShiftPrint))
            logger.debug("------------------------------------------------------------------------")
            
            TempListSequencesVector=ListResultsToListSequencesVector(NewList,EllipticConfigurationList)
            ResultList.append({'Found cnd':len(OldList[0])+1,
                               'Type':sequenceToDescriptionString(NewList[0],EllipticConfigurationList," x ", ", "),
                               'Cardinality': len(NewList),
                               'ListSequencesVector': TempListSequencesVector,
                               'ListSequencesSupport': ListVectorToListSupport(TempListSequencesVector,ShiftPrint),                         
                               'TypeVector': SequenceToVectorType(NewList[0],EllipticConfigurationList)})
            
            #Call the recursion and append the new cnd found to the list ListNewCND
            ListNewCND.append(RecursiveCnd(NewList,i,ResultList,EllipticConfigurationList,IntersectionMatrix,ShiftPrint))

        else:
            logger.debug("\t\t\t No list found. Moving on.")
            logger.debug("------------------------------------------------------------------------")
            
    return max(ListNewCND)



def SequenceToVectorType(Seq,EllipticConfigurationList):
    #SequenceToVectorType(Seq,EllipticConfigurationList) takes an isotropic sequence Seq, and returns a vector of lenght len(EllipticConfigurationList), which specifies how many elliptic fibrations of each type Seq contains
    Result=[0]*len(EllipticConfigurationList)
    for item in Seq:
        Result[item[0]]+=1
    return Result



##########################################################
# Auxiliary functions called in STEP 6
##########################################################

def IsStrictlyContained(Item,Container):
    #IsStrictlyContained(Item,Container) takes two numerical vectors of the same length, and returns true if all the entries of Item are lesser or equal than the corresponding entries of Container, AND if there is at least one entry of Item which is strictly smaller than the corresponding entry of Container. This is used on the vectors describing the type of isotropic sequences, and gives a necessary condition for a sequence of type "Item" to be contained in a sequence of type "Container"
    if (sum(Container)<=sum(Item)):
        return False
    for i in range(0,len(Item)):
        if Item[i]>Container[i]:
            return False
    return True
    
def IsSubSequence(Seq1,Seq2):
    #IsSubSequence(Seq1,Seq2) checks if all the elements of the list Seq1 are also elements of the list Seq2
    for x in Seq1:
        if not x in Seq2:
            return False
    return True

def IsInList(Seq,List):
    #IsInList(Seq,List) checks if the list Seq is contained in one of the entries of List
    for s in List:
        if IsSubSequence(Seq,s):
            return True
    return False    

    
##########################################################
#Functions to handle strings and logs
##########################################################

def ListResultsToListSequencesVector(List,Diagrams):
    #ListResultsToListSequencesVector takes in input a List of isotropic sequences, where the single elements are expressed as a pair of indices which identify the element in the list Diagrams. The output is the same list, in which the single elements are explicitly written as coefficients.
    
    ListVec=[]
    for v in List:
        TempVec=[]
        for w in v:
            TempVec.append(GetEllipticConfiguration(Diagrams,w))
        ListVec.append(TempVec)
    return ListVec

def ListVectorToListSupport(List,ShiftPrint=0):
    #ListVectorToListSupport takes a list of vectors and applies to each element the function vectorToVertexListConvert, joining all the results into a list.
    ListVec=[]
    for v in List:
        TempVec=[]
        for w in v:
            TempVec.append(vectorToVertexListConvert(w,ShiftPrint))
        ListVec.append(TempVec)
    return ListVec

def vectorToVertexListConvert(Vector,ShiftPrint=0):
    #vectorToVertexListConvert takes a vector representing an elliptic configuration as a linear combination of the original curves, and returns a list of vertices of their dual graph. The parameter "ShiftPrint" is used to specify the shift of the labels.
    #It is kind of an "inverse" of the function vertexListToVectorConvert. However, this function "forgets" the coefficients in Vector, hence its output does not distinguish between fiber and half fiber.
    
    Result=[]
    for i in range(0,len(Vector)):
        if Vector[i] != 0: Result.append(i+ShiftPrint)
    return Result

def sequenceToVerticesListToString(Seq,EllipticConfigurationList,ShiftPrint):
    #sequenceToVerticesListToString takes a list "Seq" of elliptic configurations, expressed by their indexes in EllipticConfigurationList. It first uses the function "GetEllipticConfiguration" to convert them back to vectors in terms of the original curves, and then convert each of them to a list of vertices using the function vectorToVertexListConvert, obtaining the final string. The parameter "ShiftPrint" is used to specify the shift of the labels.
    
    Result=""
    for l in Seq:
        v=GetEllipticConfiguration(EllipticConfigurationList,l)
        Result=Result+str(vectorToVertexListConvert(v,ShiftPrint))+","
    return Result[0:-1]

def sequenceToDescriptionString(Sequence,EllipticConfigurationList,Times=" x ",Plus=" + "):
    #sequenceToDescriptionString takes a isotropic sequence "Sequence" and produce a description based on the templates "N x TypeOfTheDiagram", where N is the number of elliptic configurations of type TypeOfTheDiagram appearing in the sequence. The specific names of elliptic configurations are described in "EllipticConfigurationList". It is possible to specify a custom string instead of " x " by using the optional argument "Times".
    
    if (len(Sequence)==0):
        return "[]; "
    Result=""
    Current=Sequence[0][0]
    Number=0
    for t in Sequence:
        if (t[0]==Current):
            Number=Number+1
        else:
            Result=Result + str(Number) + Times + EllipticConfigurationList[Current]['Type'] + Plus
            Current=t[0]
            Number=1
    Result=Result + str(Number) + Times + EllipticConfigurationList[Current]['Type'] + ""
    return Result

def NicePrint(FR,Basis,Matrix):
    print("--- Intersection matrix ---")
    print(Matrix.str())
    print()
    print("--- Basis of Num(S) ---")
    for i in Basis:
        print("  " + str(i))
    print()
    print("--- Value of cnd(S,R) ---")
    print("  " + str(FR['Cnd']))
    print()
    print("--- List of elliptic fibrations ---")
    for i in FR['EllipticFibrations']:
        print("  Fibrations of type " + i['Type'] + ". Total: " + str(len(i['FibrationsData'])))
        for j in i['FibrationsData']:
            print("    " + str(j['FibersSupport']))
    print()
    print("--- List of R-saturated sequences ---")
    for i in FR['SaturatedSequences']:
        print("  Sequences of type " + i['Type'])
        print("  Length " + str(i['Found cnd']))
        print("  Total number " + str(len(i['ListSequencesSupport'])))
        for j in i['ListSequencesSupport']:
            print("    " + str(j))
