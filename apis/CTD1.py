# -*- coding: utf-8 -*-
'''
#####################################################################################

This module is used for computing the composition, transition and distribution

descriptors based on the different properties of AADs. The AADs with the same

properties is marked as the same number. You can get 147 descriptors for a given

protein sequence. You can freely use and distribute it. If you hava  any problem,

you could contact with us timely!

References:

[1]: Inna Dubchak, Ilya Muchink, Stephen R.Holbrook and Sung-Hou Kim. Prediction

of protein folding class using global description of amino acid sequence. Proc.Natl.

Acad.Sci.USA, 1995, 92, 8700-8704.

[2]:Inna Dubchak, Ilya Muchink, Christopher Mayor, Igor Dralyuk and Sung-Hou Kim.

Recognition of a Protein Fold in the Context of the SCOP classification. Proteins:

Structure, Function and Genetics,1999,35,401-407.

Authors: Dongsheng Cao and Yizeng Liang.

Date: 2010.11.22

Email: oriental-cds@163.com

#####################################################################################
'''


import string, math, copy


AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

_Hydrophobicity={'1':'RKEDQN','2':'GASTPHY','3':'CLVIMFW'}
#'1'stand for Polar; '2'stand for Neutral, '3' stand for Hydrophobicity

_NormalizedVDWV={'1':'GASTPD','2':'NVEQIL','3':'MHKFRYW'}
#'1'stand for (0-2.78); '2'stand for (2.95-4.0), '3' stand for (4.03-8.08)

_Polarity={'1':'LIFWCMVY','2':'CPNVEQIL','3':'KMHFRYW'}
#'1'stand for (4.9-6.2); '2'stand for (8.0-9.2), '3' stand for (10.4-13.0)

_Charge={'1':'KR','2':'ANCQGHILMFPSTWYV','3':'DE'}
#'1'stand for Positive; '2'stand for Neutral, '3' stand for Negative

_SecondaryStr={'1':'EALMQKRH','2':'VIYCWFT','3':'GNPSD'}
#'1'stand for Helix; '2'stand for Strand, '3' stand for coil

_SolventAccessibility={'1':'ALFCGIVW','2':'RKQEND','3':'MPSTHY'}
#'1'stand for Buried; '2'stand for Exposed, '3' stand for Intermediate

_Polarizability={'1':'GASDT','2':'CPNVEQIL','3':'KMHFRYW'}
#'1'stand for (0-0.108); '2'stand for (0.128-0.186), '3' stand for (0.219-0.409)

##############################################
##         Recently added featurs           ##
##############################################

_SurfaceTension={'1':'GQDNAHR','2':'KTSEC','3':'ILMFPWYV'} # Surface tension
#'1'stand for (-0.20 ~ 0.16); '2'stand for (-0.3 ~ -0.52), '3' stand for (-0.98 ~ -2.46)

_PPIHotspotPropBogan={'1':'DHIKNPRWY','2':'EQSTGAMF','3':'CLV'} # Protein-protein interface hotspot propensity - Bogan
#'1'stand for (High (5-21%)); '2'stand for (Medium (1.12-3.64%)), '3' stand for (Low (0-0.83%))

_PPIPropMa={'1':'CDFMPQRWY','2':'AGHVLNST','3':'EIK'} # Protein-protein interface propensity - Ma
#'1'stand for (High (1.21-2.02)); '2'stand for (Medium (0.63-1.12)), '3' stand for (Low (0.14-0.29))

_PDNAIPropSchneider={'1':'GKNQRSTY','2':'ADEFHILVW','3':'CMP'} # Protein-DNA interface propensity - Schneider
#'1'stand for (High (4-30%)); '2'stand for (Medium (1-3%)), '3' stand for (Low (0-1%))

_PDNAIPropAhmad={'1':'GHKNQRSTY','2':'ADEFIPVW','3':'CLM'} # Protein-DNA interface propensity - Ahmad
#'1'stand for (High (25-100%)); '2'stand for (Medium (5-18%)), '3' stand for (Low (0-4%))

_PRNAIPropKim={'1':'HKMRY','2':'FGILNPQSVW','3':'CDEAT'} # Protein-RNA interface propensity - Kim
#'1'stand for (High (0.25-11)); '2'stand for (Medium (-0.25 –0.17)), '3' stand for (Low (-0.3 - -0.8))

_PRNAIPropEllis={'1':'HGKMRSYW','2':'AFINPQT','3':'CDELV'} # Protein-RNA interface propensity - Ellis
#'1'stand for (High (1.18-2.07)); '2'stand for (Medium (0.84-1.16)), '3' stand for (Low (0.41-0.8))

_PRNAIPropPhipps={'1':'HKMQRS','2':'ADEFGLNPVY','3':'CITW'} # Protein-RNA interface propensity - Phipps
#'1'stand for (High (0.95-1.8)); '2'stand for (Medium (0.5-0.95)), '3' stand for (Low (0-0.5))

_PLBSPropKhazanov={'1':'CFHWY','2':'GILNMSTR','3':'AEDKPQV'} # Protein-ligand binding site propensity - Khazanov
#'1'stand for (High (≥2.25)); '2'stand for (Medium (1.6-2.3)), '3' stand for (Low (≤1.5))

_PLVBSKhazanov={'1':'CFHWYM','2':'DGILNSTV','3':'AEKPQR'} # Protein-ligand valid binding site propensity - Khazanov
#'1'stand for (High (≥1.4)); '2'stand for (Medium (0.79-1.21)), '3' stand for (Low (≤0.76))

_PropPLPANBIntImai={'1':'DEHRY','2':'CFKMNQSTW','3':'AGILPV'} # Propensity for protein-ligand polar & aromatic non-bonded interactions-Imai
#'1'stand for (High (477-1197)); '2'stand for (Medium (95-423)), '3' stand for (Low (<95))

_MolecularWeight={'1':'AGS','2':'CDEHIKLMNQPTV','3':'FRWY'} # Molecular Weight
#'1'stand for (Low (75-105)); '2'stand for (Medium(115-155)), '3' stand for (High (165-204))

_cLogP={'1':'RKDNEQH','2':'PYSTGACV','3':'WMFLI'} # cLogP
#'1'stand for (-4.2 - -3.3); '2'stand for (-3.07 – 2.26), '3' stand for (-1.78 - -1.05)

_NoHydroBondDonorSideChain={'1':'HKNQR','2':'DESTWY','3':'ACGFILMPV'} # No of hydrogen bond donor in side chain
#'1'stand for (>1); '2'stand for (1), '3' stand for (0)

_NoHydroBondAccSideChain={'1':'DEHNQR','2':'KSTWY','3':'ACGFILMPV'} # No of hydrogen bond acceptor in side chain
#'1'stand for (>1); '2'stand for (1), '3' stand for (0)

_SolubilityInWater={'1':'ACGKRT','2':'EFHILMNPQSVW','3':'DY'} # Solubility in water
#'1'stand for (High (9-65 g/100g)); '2'stand for (Medium (1.14-7.44 g/100g)), '3' stand for (Low (0.048-0.82g/100g))

_AminoAcidFlexInd={'1':'EGKNQS','2':'ADHIPRTV','3':'CFLMWY'} # Amino acid flexibility index
#'1'stand for (Very flexible); '2'stand for (Moderately flexible), '3' stand for (Less flexible)


##You can continuely add other properties of AADs to compute descriptors of protein sequence.

_AATProperty=(_Hydrophobicity,_NormalizedVDWV,_Polarity,_Charge,_SecondaryStr,_SolventAccessibility,_Polarizability,\
_SurfaceTension,_PPIHotspotPropBogan,_PPIPropMa,_PDNAIPropSchneider,_PDNAIPropAhmad,_PRNAIPropKim,_PRNAIPropEllis,_PRNAIPropPhipps,\
_PRNAIPropPhipps,_PLBSPropKhazanov,_PLVBSKhazanov,_PropPLPANBIntImai,_MolecularWeight,_cLogP,_NoHydroBondDonorSideChain,\
_NoHydroBondAccSideChain,_SolubilityInWater,_AminoAcidFlexInd)

_AATPropertyName=('_Hydrophobicity','_NormalizedVDWV','_Polarity','_Charge','_SecondaryStr','_SolventAccessibility','_Polarizability',\
'_SurfaceTension','_PPIHotspotPropBogan','_PPIPropMa','_PDNAIPropSchneider','_PDNAIPropAhmad','_PRNAIPropKim','_PRNAIPropEllis','_PRNAIPropPhipps',\
'_PRNAIPropPhipps','_PLBSPropKhazanov','_PLVBSKhazanov','_PropPLPANBIntImai','_MolecularWeight','_cLogP','_NoHydroBondDonorSideChain',\
'_NoHydroBondAccSideChain','_SolubilityInWater','_AminoAcidFlexInd')


##################################################################################################

def StringtoNum(ProteinSequence,AAProperty):
	"""
	###############################################################################################
	Tranform the protein sequence into the string form such as 32123223132121123.

	Usage:

	result=StringtoNum(protein,AAProperty)

	Input: protein is a pure protein sequence.

	AAProperty is a dict form containing classifciation of amino acids such as _Polarizability.

	Output: result is a string such as 123321222132111123222
	###############################################################################################
	"""

	hardProteinSequence=copy.deepcopy(ProteinSequence)
	for k,m in AAProperty.items():
		for index in m:
			hardProteinSequence=str.replace(hardProteinSequence,index,k)
	TProteinSequence=hardProteinSequence

	return TProteinSequence


def CalculateComposition(ProteinSequence,AAProperty,AAPName):
	"""
	###############################################################################################
	A method used for computing composition descriptors.

	Usage:

	result=CalculateComposition(protein,AAProperty,AAPName)

	Input: protein is a pure protein sequence.

	AAProperty is a dict form containing classifciation of amino acids such as _Polarizability.

	AAPName is a string used for indicating a AAP name.

	Output: result is a dict form containing composition descriptors based on the given property.
	###############################################################################################
	"""
	TProteinSequence=StringtoNum(ProteinSequence,AAProperty)
	Result={}
	Num=len(TProteinSequence)
	Result[AAPName+'C'+'1']=round(float(TProteinSequence.count('1'))/Num,3)
	Result[AAPName+'C'+'2']=round(float(TProteinSequence.count('2'))/Num,3)
	Result[AAPName+'C'+'3']=round(float(TProteinSequence.count('3'))/Num,3)
	return Result

def CalculateTransition(ProteinSequence,AAProperty,AAPName):
	"""
	###############################################################################################
	A method used for computing transition descriptors

	Usage:

	result=CalculateTransition(protein,AAProperty,AAPName)

	Input:protein is a pure protein sequence.

	AAProperty is a dict form containing classifciation of amino acids such as _Polarizability.

	AAPName is a string used for indicating a AAP name.

	Output:result is a dict form containing transition descriptors based on the given property.
	###############################################################################################
	"""

	TProteinSequence=StringtoNum(ProteinSequence,AAProperty)
	Result={}
	Num=len(TProteinSequence)
	CTD=TProteinSequence
	Result[AAPName+'T'+'12']=round(float(CTD.count('12')+CTD.count('21'))/(Num-1),3)
	Result[AAPName+'T'+'13']=round(float(CTD.count('13')+CTD.count('31'))/(Num-1),3)
	Result[AAPName+'T'+'23']=round(float(CTD.count('23')+CTD.count('32'))/(Num-1),3)
	return Result

def CalculateDistribution(ProteinSequence,AAProperty,AAPName):

	"""
	###############################################################################################
	A method used for computing distribution descriptors.

	Usage:

	result=CalculateDistribution(protein,AAProperty,AAPName)

	Input:protein is a pure protein sequence.

	AAProperty is a dict form containing classifciation of amino acids such as _Polarizability.

	AAPName is a string used for indicating a AAP name.

	Output:result is a dict form containing Distribution descriptors based on the given property.
	###############################################################################################
	"""
	TProteinSequence=StringtoNum(ProteinSequence,AAProperty)
	Result={}
	Num=len(TProteinSequence)
	temp=('1','2','3')
	for i in temp:
		num=TProteinSequence.count(i)
		ink=1
		indexk=0
		cds=[]
		while ink<=num:
			indexk=str.find(TProteinSequence,i,indexk)+1
			cds.append(indexk)
			ink=ink+1

		if cds==[]:
			Result[AAPName+'D'+i+'001']=0
			Result[AAPName+'D'+i+'025']=0
			Result[AAPName+'D'+i+'050']=0
			Result[AAPName+'D'+i+'075']=0
			Result[AAPName+'D'+i+'100']=0
		else:

			Result[AAPName+'D'+i+'001']=round(float(cds[0])/Num*100,3)
			Result[AAPName+'D'+i+'025']=round(float(cds[int(math.floor(num*0.25))-1])/Num*100,3)
			Result[AAPName+'D'+i+'050']=round(float(cds[int(math.floor(num*0.5))-1])/Num*100,3)
			Result[AAPName+'D'+i+'075']=round(float(cds[int(math.floor(num*0.75))-1])/Num*100,3)
			Result[AAPName+'D'+i+'100']=round(float(cds[-1])/Num*100,3)

	return Result

##################################################################################################
def CalculateCompositionHydrophobicity(ProteinSequence):

	"""
	###############################################################################################
	A method used for calculating composition descriptors based on Hydrophobicity of

	AADs.

	Usage:

	result=CalculateCompositionHydrophobicity(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Composition descriptors based on Hydrophobicity.
	###############################################################################################
	"""

	result=CalculateComposition(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
	return result

def CalculateCompositionNormalizedVDWV(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating composition descriptors based on NormalizedVDWV of

	AADs.

	Usage:

	result=CalculateCompositionNormalizedVDWV(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Composition descriptors based on NormalizedVDWV.
	###############################################################################################
	"""
	result=CalculateComposition(ProteinSequence,_NormalizedVDWV,'_NormalizedVDWV')
	return result

def CalculateCompositionPolarity(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating composition descriptors based on Polarity of

	AADs.

	Usage:

	result=CalculateCompositionPolarity(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Composition descriptors based on Polarity.
	###############################################################################################
	"""

	result=CalculateComposition(ProteinSequence,_Polarity,'_Polarity')
	return result

def CalculateCompositionCharge(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating composition descriptors based on Charge of

	AADs.

	Usage:

	result=CalculateCompositionCharge(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Composition descriptors based on Charge.
	###############################################################################################
	"""

	result=CalculateComposition(ProteinSequence,_Charge,'_Charge')
	return result

def CalculateCompositionSecondaryStr(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating composition descriptors based on SecondaryStr of

	AADs.

	Usage:

	result=CalculateCompositionSecondaryStr(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Composition descriptors based on SecondaryStr.
	###############################################################################################
	"""

	result=CalculateComposition(ProteinSequence,_SecondaryStr,'_SecondaryStr')
	return result

def CalculateCompositionSolventAccessibility(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating composition descriptors based on SolventAccessibility

	of  AADs.

	Usage:

	result=CalculateCompositionSolventAccessibility(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Composition descriptors based on SolventAccessibility.
	###############################################################################################
	"""

	result=CalculateComposition(ProteinSequence,_SolventAccessibility,'_SolventAccessibility')
	return result

def CalculateCompositionPolarizability(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating composition descriptors based on Polarizability of

	AADs.

	Usage:

	result=CalculateCompositionPolarizability(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Composition descriptors based on Polarizability.
	###############################################################################################
	"""

	result=CalculateComposition(ProteinSequence,_Polarizability,'_Polarizability')
	return result

def CalculateCompositionSurfaceTension(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_SurfaceTension,'_SurfaceTension')
	return result

def CalculateCompositionPPIHotspotPropBogan(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_PPIHotspotPropBogan,'_PPIHotspotPropBogan')
	return result

def CalculateCompositionPPIPropMa(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_PPIPropMa,'_PPIPropMa')
	return result

def CalculateCompositionPDNAIPropSchneider(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_PDNAIPropSchneider,'_PDNAIPropSchneider')
	return result

def CalculateCompositionPDNAIPropAhmad(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_PDNAIPropAhmad,'_PDNAIPropAhmad')
	return result

def CalculateCompositionPRNAIPropKim(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_PRNAIPropKim,'_PRNAIPropKim')
	return result

def CalculateCompositionPRNAIPropEllis(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_PRNAIPropEllis,'_PRNAIPropEllis')
	return result

def CalculateCompositionPRNAIPropPhipps(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_PRNAIPropPhipps,'_PRNAIPropPhipps')
	return result

def CalculateCompositionPRNAIPropPhipps(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_PRNAIPropPhipps,'_PRNAIPropPhipps')
	return result

def CalculateCompositionPLBSPropKhazanov(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_PLBSPropKhazanov,'_PLBSPropKhazanov')
	return result

def CalculateCompositionPLVBSKhazanov(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_PLVBSKhazanov,'_PLVBSKhazanov')
	return result

def CalculateCompositionPropPLPANBIntImai(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_PropPLPANBIntImai,'_PropPLPANBIntImai')
	return result

def CalculateCompositionMolecularWeight(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_MolecularWeight,'_MolecularWeight')
	return result

def CalculateCompositioncLogP(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_cLogP,'_cLogP')
	return result

def CalculateCompositionNoHydroBondDonorSideChain(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_NoHydroBondDonorSideChain,'_NoHydroBondDonorSideChain')
	return result

def CalculateCompositionNoHydroBondAccSideChain(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_NoHydroBondAccSideChain,'_NoHydroBondAccSideChain')
	return result

def CalculateCompositionSolubilityInWater(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_SolubilityInWater,'_SolubilityInWater')
	return result

def CalculateCompositionAminoAcidFlexInd(ProteinSequence):
	result=CalculateComposition(ProteinSequence,_AminoAcidFlexInd,'_AminoAcidFlexInd')
	return result

##################################################################################################


##################################################################################################
def CalculateTransitionHydrophobicity(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating Transition descriptors based on Hydrophobicity of

	AADs.

	Usage:

	result=CalculateTransitionHydrophobicity(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Transition descriptors based on Hydrophobicity.
	###############################################################################################
	"""

	result=CalculateTransition(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
	return result

def CalculateTransitionNormalizedVDWV(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating Transition descriptors based on NormalizedVDWV of

	AADs.

	Usage:

	result=CalculateTransitionNormalizedVDWV(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Transition descriptors based on NormalizedVDWV.
	###############################################################################################
	"""

	result=CalculateTransition(ProteinSequence,_NormalizedVDWV,'_NormalizedVDWV')
	return result

def CalculateTransitionPolarity(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating Transition descriptors based on Polarity of

	AADs.

	Usage:

	result=CalculateTransitionPolarity(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Transition descriptors based on Polarity.
	###############################################################################################
	"""

	result=CalculateTransition(ProteinSequence,_Polarity,'_Polarity')
	return result

def CalculateTransitionCharge(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating Transition descriptors based on Charge of

	AADs.

	Usage:

	result=CalculateTransitionCharge(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Transition descriptors based on Charge.
	###############################################################################################
	"""

	result=CalculateTransition(ProteinSequence,_Charge,'_Charge')
	return result

def CalculateTransitionSecondaryStr(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating Transition descriptors based on SecondaryStr of

	AADs.

	Usage:

	result=CalculateTransitionSecondaryStr(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Transition descriptors based on SecondaryStr.
	###############################################################################################
	"""

	result=CalculateTransition(ProteinSequence,_SecondaryStr,'_SecondaryStr')
	return result

def CalculateTransitionSolventAccessibility(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating Transition descriptors based on SolventAccessibility

	of  AADs.

	Usage:

	result=CalculateTransitionSolventAccessibility(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Transition descriptors based on SolventAccessibility.
	###############################################################################################
	"""

	result=CalculateTransition(ProteinSequence,_SolventAccessibility,'_SolventAccessibility')
	return result

def CalculateTransitionPolarizability(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating Transition descriptors based on Polarizability of

	AADs.

	Usage:

	result=CalculateTransitionPolarizability(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Transition descriptors based on Polarizability.
	###############################################################################################
	"""

	result=CalculateTransition(ProteinSequence,_Polarizability,'_Polarizability')
	return result

def CalculateTransitionSurfaceTension(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_SurfaceTension,'_SurfaceTension')
	return result

def CalculateTransitionPPIHotspotPropBogan(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_PPIHotspotPropBogan,'_PPIHotspotPropBogan')
	return result

def CalculateTransitionPPIPropMa(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_PPIPropMa,'_PPIPropMa')
	return result

def CalculateTransitionPDNAIPropSchneider(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_PDNAIPropSchneider,'_PDNAIPropSchneider')
	return result

def CalculateTransitionPDNAIPropAhmad(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_PDNAIPropAhmad,'_PDNAIPropAhmad')
	return result

def CalculateTransitionPRNAIPropKim(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_PRNAIPropKim,'_PRNAIPropKim')
	return result

def CalculateTransitionPRNAIPropEllis(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_PRNAIPropEllis,'_PRNAIPropEllis')
	return result

def CalculateTransitionPRNAIPropPhipps(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_PRNAIPropPhipps,'_PRNAIPropPhipps')
	return result

def CalculateTransitionPRNAIPropPhipps(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_PRNAIPropPhipps,'_PRNAIPropPhipps')
	return result

def CalculateTransitionPLBSPropKhazanov(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_PLBSPropKhazanov,'_PLBSPropKhazanov')
	return result

def CalculateTransitionPLVBSKhazanov(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_PLVBSKhazanov,'_PLVBSKhazanov')
	return result

def CalculateTransitionPropPLPANBIntImai(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_PropPLPANBIntImai,'_PropPLPANBIntImai')
	return result

def CalculateTransitionMolecularWeight(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_MolecularWeight,'_MolecularWeight')
	return result

def CalculateTransitioncLogP(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_cLogP,'_cLogP')
	return result

def CalculateTransitionNoHydroBondDonorSideChain(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_NoHydroBondDonorSideChain,'_NoHydroBondDonorSideChain')
	return result

def CalculateTransitionNoHydroBondAccSideChain(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_NoHydroBondAccSideChain,'_NoHydroBondAccSideChain')
	return result

def CalculateTransitionSolubilityInWater(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_SolubilityInWater,'_SolubilityInWater')
	return result

def CalculateTransitionAminoAcidFlexInd(ProteinSequence):
	result=CalculateTransition(ProteinSequence,_AminoAcidFlexInd,'_AminoAcidFlexInd')
	return result

##################################################################################################
##################################################################################################
def CalculateDistributionHydrophobicity(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating Distribution descriptors based on Hydrophobicity of

	AADs.

	Usage:

	result=CalculateDistributionHydrophobicity(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Distribution descriptors based on Hydrophobicity.
	###############################################################################################
	"""

	result=CalculateDistribution(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
	return result

def CalculateDistributionNormalizedVDWV(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating Distribution descriptors based on NormalizedVDWV of

	AADs.

	Usage:

	result=CalculateDistributionNormalizedVDWV(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Distribution descriptors based on NormalizedVDWV.
	###############################################################################################
	"""

	result=CalculateDistribution(ProteinSequence,_NormalizedVDWV,'_NormalizedVDWV')
	return result

def CalculateDistributionPolarity(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating Distribution descriptors based on Polarity of

	AADs.

	Usage:

	result=CalculateDistributionPolarity(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Distribution descriptors based on Polarity.
	###############################################################################################
	"""

	result=CalculateDistribution(ProteinSequence,_Polarity,'_Polarity')
	return result

def CalculateDistributionCharge(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating Distribution descriptors based on Charge of

	AADs.

	Usage:

	result=CalculateDistributionCharge(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Distribution descriptors based on Charge.
	###############################################################################################
	"""

	result=CalculateDistribution(ProteinSequence,_Charge,'_Charge')
	return result

def CalculateDistributionSecondaryStr(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating Distribution descriptors based on SecondaryStr of

	AADs.

	Usage:

	result=CalculateDistributionSecondaryStr(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Distribution descriptors based on SecondaryStr.
	###############################################################################################
	"""

	result=CalculateDistribution(ProteinSequence,_SecondaryStr,'_SecondaryStr')
	return result

def CalculateDistributionSolventAccessibility(ProteinSequence):

	"""
	###############################################################################################
	A method used for calculating Distribution descriptors based on SolventAccessibility

	of  AADs.

	Usage:

	result=CalculateDistributionSolventAccessibility(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Distribution descriptors based on SolventAccessibility.
	###############################################################################################
	"""

	result=CalculateDistribution(ProteinSequence,_SolventAccessibility,'_SolventAccessibility')
	return result

def CalculateDistributionPolarizability(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating Distribution descriptors based on Polarizability of

	AADs.

	Usage:

	result=CalculateDistributionPolarizability(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing Distribution descriptors based on Polarizability.
	###############################################################################################
	"""

	result=CalculateDistribution(ProteinSequence,_Polarizability,'_Polarizability')
	return result

def CalculateDistributionSurfaceTension(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_SurfaceTension,'_SurfaceTension')
	return result

def CalculateDistributionPPIHotspotPropBogan(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_PPIHotspotPropBogan,'_PPIHotspotPropBogan')
	return result

def CalculateDistributionPPIPropMa(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_PPIPropMa,'_PPIPropMa')
	return result

def CalculateDistributionPDNAIPropSchneider(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_PDNAIPropSchneider,'_PDNAIPropSchneider')
	return result

def CalculateDistributionPDNAIPropAhmad(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_PDNAIPropAhmad,'_PDNAIPropAhmad')
	return result

def CalculateDistributionPRNAIPropKim(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_PRNAIPropKim,'_PRNAIPropKim')
	return result

def CalculateDistributionPRNAIPropEllis(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_PRNAIPropEllis,'_PRNAIPropEllis')
	return result

def CalculateDistributionPRNAIPropPhipps(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_PRNAIPropPhipps,'_PRNAIPropPhipps')
	return result

def CalculateDistributionPRNAIPropPhipps(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_PRNAIPropPhipps,'_PRNAIPropPhipps')
	return result

def CalculateDistributionPLBSPropKhazanov(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_PLBSPropKhazanov,'_PLBSPropKhazanov')
	return result

def CalculateDistributionPLVBSKhazanov(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_PLVBSKhazanov,'_PLVBSKhazanov')
	return result

def CalculateDistributionPropPLPANBIntImai(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_PropPLPANBIntImai,'_PropPLPANBIntImai')
	return result

def CalculateDistributionMolecularWeight(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_MolecularWeight,'_MolecularWeight')
	return result

def CalculateDistributioncLogP(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_cLogP,'_cLogP')
	return result

def CalculateDistributionNoHydroBondDonorSideChain(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_NoHydroBondDonorSideChain,'_NoHydroBondDonorSideChain')
	return result

def CalculateDistributionNoHydroBondAccSideChain(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_NoHydroBondAccSideChain,'_NoHydroBondAccSideChain')
	return result

def CalculateDistributionSolubilityInWater(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_SolubilityInWater,'_SolubilityInWater')
	return result

def CalculateDistributionAminoAcidFlexInd(ProteinSequence):
	result=CalculateDistribution(ProteinSequence,_AminoAcidFlexInd,'_AminoAcidFlexInd')
	return result

##################################################################################################

def CalculateC(ProteinSequence):
	"""
	###############################################################################################
	Calculate all composition descriptors based seven different properties of AADs.

	Usage:

	result=CalculateC(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing all composition descriptors.
	###############################################################################################
	"""
	result={}
	result.update(CalculateCompositionPolarizability(ProteinSequence))
	result.update(CalculateCompositionSolventAccessibility(ProteinSequence))
	result.update(CalculateCompositionSecondaryStr(ProteinSequence))
	result.update(CalculateCompositionCharge(ProteinSequence))
	result.update(CalculateCompositionPolarity(ProteinSequence))
	result.update(CalculateCompositionNormalizedVDWV(ProteinSequence))
	result.update(CalculateCompositionHydrophobicity(ProteinSequence))
	return result

def CalculateT(ProteinSequence):
	"""
	###############################################################################################
	Calculate all transition descriptors based seven different properties of AADs.

	Usage:

	result=CalculateT(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing all transition descriptors.
	###############################################################################################
	"""
	result={}
	result.update(CalculateTransitionPolarizability(ProteinSequence))
	result.update(CalculateTransitionSolventAccessibility(ProteinSequence))
	result.update(CalculateTransitionSecondaryStr(ProteinSequence))
	result.update(CalculateTransitionCharge(ProteinSequence))
	result.update(CalculateTransitionPolarity(ProteinSequence))
	result.update(CalculateTransitionNormalizedVDWV(ProteinSequence))
	result.update(CalculateTransitionHydrophobicity(ProteinSequence))
	return result

def CalculateD(ProteinSequence):
	"""
	###############################################################################################
	Calculate all distribution descriptors based seven different properties of AADs.

	Usage:

	result=CalculateD(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing all distribution descriptors.
	###############################################################################################
	"""
	result={}
	result.update(CalculateDistributionPolarizability(ProteinSequence))
	result.update(CalculateDistributionSolventAccessibility(ProteinSequence))
	result.update(CalculateDistributionSecondaryStr(ProteinSequence))
	result.update(CalculateDistributionCharge(ProteinSequence))
	result.update(CalculateDistributionPolarity(ProteinSequence))
	result.update(CalculateDistributionNormalizedVDWV(ProteinSequence))
	result.update(CalculateDistributionHydrophobicity(ProteinSequence))
	return result


def CalculateCTD(ProteinSequence, output_file_path=''):
	"""
	###############################################################################################
	Calculate all CTD descriptors based seven different properties of AADs.

	Usage:

	result=CalculateCTD(protein)

	Input:protein is a pure protein sequence.

	Output:result is a dict form containing all CTD descriptors.
	###############################################################################################
	"""
	result={}
	## Calculate Composition
	result.update(CalculateCompositionPolarizability(ProteinSequence))
	result.update(CalculateCompositionSolventAccessibility(ProteinSequence))
	result.update(CalculateCompositionSecondaryStr(ProteinSequence))
	result.update(CalculateCompositionCharge(ProteinSequence))
	result.update(CalculateCompositionPolarity(ProteinSequence))
	result.update(CalculateCompositionNormalizedVDWV(ProteinSequence))
	result.update(CalculateCompositionHydrophobicity(ProteinSequence))
	result.update(CalculateCompositionSurfaceTension(ProteinSequence))
	result.update(CalculateCompositionPPIHotspotPropBogan(ProteinSequence))
	result.update(CalculateCompositionPPIPropMa(ProteinSequence))
	result.update(CalculateCompositionPDNAIPropSchneider(ProteinSequence))
	result.update(CalculateCompositionPDNAIPropAhmad(ProteinSequence))
	result.update(CalculateCompositionPRNAIPropKim(ProteinSequence))
	result.update(CalculateCompositionPRNAIPropEllis(ProteinSequence))
	result.update(CalculateCompositionPRNAIPropPhipps(ProteinSequence))
	result.update(CalculateCompositionPLBSPropKhazanov(ProteinSequence))
	result.update(CalculateCompositionPLVBSKhazanov(ProteinSequence))
	result.update(CalculateCompositionPropPLPANBIntImai(ProteinSequence))
	result.update(CalculateCompositionMolecularWeight(ProteinSequence))
	result.update(CalculateCompositioncLogP(ProteinSequence))
	result.update(CalculateCompositionNoHydroBondDonorSideChain(ProteinSequence))
	result.update(CalculateCompositionNoHydroBondAccSideChain(ProteinSequence))
	result.update(CalculateCompositionSolubilityInWater(ProteinSequence))
	result.update(CalculateCompositionAminoAcidFlexInd(ProteinSequence))

	## Calculate Transition
	result.update(CalculateTransitionPolarizability(ProteinSequence))
	result.update(CalculateTransitionSolventAccessibility(ProteinSequence))
	result.update(CalculateTransitionSecondaryStr(ProteinSequence))
	result.update(CalculateTransitionCharge(ProteinSequence))
	result.update(CalculateTransitionPolarity(ProteinSequence))
	result.update(CalculateTransitionNormalizedVDWV(ProteinSequence))
	result.update(CalculateTransitionHydrophobicity(ProteinSequence))
	result.update(CalculateTransitionSurfaceTension(ProteinSequence))
	result.update(CalculateTransitionPPIHotspotPropBogan(ProteinSequence))
	result.update(CalculateTransitionPPIPropMa(ProteinSequence))
	result.update(CalculateTransitionPDNAIPropSchneider(ProteinSequence))
	result.update(CalculateTransitionPDNAIPropAhmad(ProteinSequence))
	result.update(CalculateTransitionPRNAIPropKim(ProteinSequence))
	result.update(CalculateTransitionPRNAIPropEllis(ProteinSequence))
	result.update(CalculateTransitionPRNAIPropPhipps(ProteinSequence))
	result.update(CalculateTransitionPLBSPropKhazanov(ProteinSequence))
	result.update(CalculateTransitionPLVBSKhazanov(ProteinSequence))
	result.update(CalculateTransitionPropPLPANBIntImai(ProteinSequence))
	result.update(CalculateTransitionMolecularWeight(ProteinSequence))
	result.update(CalculateTransitioncLogP(ProteinSequence))
	result.update(CalculateTransitionNoHydroBondDonorSideChain(ProteinSequence))
	result.update(CalculateTransitionNoHydroBondAccSideChain(ProteinSequence))
	result.update(CalculateTransitionSolubilityInWater(ProteinSequence))
	result.update(CalculateTransitionAminoAcidFlexInd(ProteinSequence))

	## Calculate Distribution
	result.update(CalculateDistributionPolarizability(ProteinSequence))
	result.update(CalculateDistributionSolventAccessibility(ProteinSequence))
	result.update(CalculateDistributionSecondaryStr(ProteinSequence))
	result.update(CalculateDistributionCharge(ProteinSequence))
	result.update(CalculateDistributionPolarity(ProteinSequence))
	result.update(CalculateDistributionNormalizedVDWV(ProteinSequence))
	result.update(CalculateDistributionHydrophobicity(ProteinSequence))
	result.update(CalculateDistributionSurfaceTension(ProteinSequence))
	result.update(CalculateDistributionPPIHotspotPropBogan(ProteinSequence))
	result.update(CalculateDistributionPPIPropMa(ProteinSequence))
	result.update(CalculateDistributionPDNAIPropSchneider(ProteinSequence))
	result.update(CalculateDistributionPDNAIPropAhmad(ProteinSequence))
	result.update(CalculateDistributionPRNAIPropKim(ProteinSequence))
	result.update(CalculateDistributionPRNAIPropEllis(ProteinSequence))
	result.update(CalculateDistributionPRNAIPropPhipps(ProteinSequence))
	result.update(CalculateDistributionPLBSPropKhazanov(ProteinSequence))
	result.update(CalculateDistributionPLVBSKhazanov(ProteinSequence))
	result.update(CalculateDistributionPropPLPANBIntImai(ProteinSequence))
	result.update(CalculateDistributionMolecularWeight(ProteinSequence))
	result.update(CalculateDistributioncLogP(ProteinSequence))
	result.update(CalculateDistributionNoHydroBondDonorSideChain(ProteinSequence))
	result.update(CalculateDistributionNoHydroBondAccSideChain(ProteinSequence))
	result.update(CalculateDistributionSolubilityInWater(ProteinSequence))
	result.update(CalculateDistributionAminoAcidFlexInd(ProteinSequence))


	if len(output_file_path)>0:
		with open(output_file_path, 'w') as f:
			f.write(str(result))

	return result
##################################################################################################

def CalculateCTD4All(input_file_path, output_file_path = ""):
	with open(output_file_path, 'w') as f:
		pass

	lines = []
	with open(input_file_path) as f:
		lines = f.readlines()

	keys = list(CalculateCTD("AD"))
	# with open(output_file_path, 'a') as f:
	# 	keys = list(CalculateCTD("AD"))
	# 	f.write("seq,")
	# 	f.write(str(keys).replace("[","").replace("]","") + ",\n")

	# print(lines[0])
	output_results = []
	for pep in lines:
		caledCTD = CalculateCTD(pep)
		results = str(pep.split('\n')[0] + ",")
		for key in keys:
			results = results + str(caledCTD.get(key))+","
		results = results + "\n"
		output_results.append(results)

	if output_file_path != "":
		with open(output_file_path, 'w') as f:
			f.write("seq,")
			f.write(str(keys).replace("[","").replace("]","").replace("'","").replace("_","") + ",\n")
			f.writelines(output_results)
	else:
		return output_results




if __name__=="__main__":

#	import scipy,string

#	result=scipy.zeros((268,147))
#	f=file('protein1.txt','r')
#	for i,j in enumerate(f:
#		temp=CalculateCTD(string.strip(j))
#		result[i,:]=temp.values()
#	scipy.savetxt('ResultNCTRER.txt', result, fmt='%15.5f',delimiter='')
#
	protein="ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"
#	print StringtoNum(protein,_Hydrophobicity)
#	print CalculateComposition(protein,_Hydrophobicity,'_Hydrophobicity')
#	print CalculateTransition(protein,_Hydrophobicity,'_Hydrophobicity')
#	print CalculateDistribution(protein,_Hydrophobicity,'_Hydrophobicity')
#	print CalculateDistributionSolventAccessibility(protein)
#	print len(CalculateCTD(protein))
#	print len(CalculateC(protein))
#	print len(CalculateT(protein))
#	print len(CalculateD(protein))
	print (CalculateCTD(protein,''))
