import kinetics
import matplotlib.pyplot as plt
from scipy.stats import norm

class vPts1_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vPts1"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kPts1','EIIAtotal','kmPts1']
		self.reaction_substrate_names = ['cAMP','PEP','PYR','pH','EIIAP','FBP','ATP']
		self.substrates = ['PEP']
		self.products = ['PYR','EIIAP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		PEP = substrates[1]
		PYR = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		FBP = substrates[5]
		ATP = substrates[6]


		kPts1 = parameters[0]
		EIIAtotal = parameters[1]
		kmPts1 = parameters[2]



		EIIA=EIIAtotal-EIIAP
		rate=kPts1*PEP*EIIA-kmPts1*PYR*EIIAP
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vPts4_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vPts4"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['vPts4_max','KPts_EIIA','KPts_GLC']
		self.reaction_substrate_names = ['EIIAP','GLCex']
		self.substrates = ['EIIAP']
		self.products = ['G6P']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		EIIAP = substrates[0]
		GLCex = substrates[1]


		vPts4_max = parameters[0]
		KPts_EIIA = parameters[1]
		KPts_GLC = parameters[2]



		rate=vPts4_max*EIIAP*GLCex/((KPts_EIIA+EIIAP)*(KPts_GLC+GLCex))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vPts4_medium_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vPts4_medium"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['vPts4_max','KPts_EIIA','rho','KPts_GLC']
		self.reaction_substrate_names = ['EIIAP','X','GLCex']
		self.substrates = ['GLCex']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		EIIAP = substrates[0]
		X = substrates[1]
		GLCex = substrates[2]


		vPts4_max = parameters[0]
		KPts_EIIA = parameters[1]
		rho = parameters[2]
		KPts_GLC = parameters[3]



		vPts4=vPts4_max*EIIAP*GLCex/((KPts_EIIA+EIIAP)*(KPts_GLC+GLCex))
		rate=vPts4*X/rho
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vNonpts_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vNonpts"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KNonpts_I','KNonpts_S','vNonpts_max','EIIAtotal']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','PYR','ATP','GLCex']
		self.substrates = []
		self.products = ['GLC']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		PYR = substrates[4]
		ATP = substrates[5]
		GLCex = substrates[6]


		KNonpts_I = parameters[0]
		KNonpts_S = parameters[1]
		vNonpts_max = parameters[2]
		EIIAtotal = parameters[3]



		EIIA=EIIAtotal-EIIAP
		rate=vNonpts_max*GLCex/(KNonpts_S+(1+EIIA/KNonpts_I)*GLCex)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vNonpts_medium_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vNonpts_medium"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KNonpts_I','KNonpts_S','rho','vNonpts_max','EIIAtotal']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','PYR','X','ATP','GLCex']
		self.substrates = ['GLCex']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		PYR = substrates[4]
		X = substrates[5]
		ATP = substrates[6]
		GLCex = substrates[7]


		KNonpts_I = parameters[0]
		KNonpts_S = parameters[1]
		rho = parameters[2]
		vNonpts_max = parameters[3]
		EIIAtotal = parameters[4]



		EIIA=EIIAtotal-EIIAP
		vNonpts=vNonpts_max*GLCex/(KNonpts_S+(1+EIIA/KNonpts_I)*GLCex)
		rate=vNonpts*X/rho
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Glk_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Glk"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KGlk_G6P_i','KGlk_GLC_m','kGlk_cat','KGlk_ATP_m']
		self.reaction_substrate_names = ['GLC','G6P','ATP','Glk']
		self.substrates = ['GLC','ATP']
		self.products = ['G6P','ADP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		GLC = substrates[0]
		G6P = substrates[1]
		ATP = substrates[2]
		Glk = substrates[3]


		KGlk_G6P_i = parameters[0]
		KGlk_GLC_m = parameters[1]
		kGlk_cat = parameters[2]
		KGlk_ATP_m = parameters[3]



		rate=Glk*kGlk_cat*(GLC/KGlk_GLC_m)*(ATP/(KGlk_ATP_m*(1+G6P/KGlk_G6P_i)))/(1+GLC/KGlk_GLC_m+ATP/(KGlk_ATP_m*(1+G6P/KGlk_G6P_i))+GLC*ATP/(KGlk_GLC_m*KGlk_ATP_m*(1+G6P/KGlk_G6P_i))+G6P/KGlk_G6P_i)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Pgi_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Pgi"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KPgi_eq','KPgi_F6P','KPgi_F6P_6pginh','KPgi_G6P_6pginh','KPgi_G6P','vPgi_max']
		self.reaction_substrate_names = ['G6P','F6P','sixPG']
		self.substrates = ['G6P']
		self.products = ['F6P']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		G6P = substrates[0]
		F6P = substrates[1]
		sixPG = substrates[2]


		KPgi_eq = parameters[0]
		KPgi_F6P = parameters[1]
		KPgi_F6P_6pginh = parameters[2]
		KPgi_G6P_6pginh = parameters[3]
		KPgi_G6P = parameters[4]
		vPgi_max = parameters[5]



		rate=vPgi_max*(G6P-F6P/KPgi_eq)/(KPgi_G6P*(1+F6P/(KPgi_F6P*(1+sixPG/KPgi_F6P_6pginh))+sixPG/KPgi_G6P_6pginh)+G6P)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Pfk_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Pfk"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KPfk_ADP_a','kPfk_cat','KPfk_AMP_b','KPfk_ADP_b','KPfk_F6P_s','KPfk_ATP_s','KPfk_PEP','LPfk','KPfk_ADP_c','KPfk_AMP_a','nPfk']
		self.reaction_substrate_names = ['F6P','PEP','Pfk','ATP','ADP','AMP']
		self.substrates = ['F6P','ATP']
		self.products = ['FBP','ADP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		F6P = substrates[0]
		PEP = substrates[1]
		Pfk = substrates[2]
		ATP = substrates[3]
		ADP = substrates[4]
		AMP = substrates[5]


		KPfk_ADP_a = parameters[0]
		kPfk_cat = parameters[1]
		KPfk_AMP_b = parameters[2]
		KPfk_ADP_b = parameters[3]
		KPfk_F6P_s = parameters[4]
		KPfk_ATP_s = parameters[5]
		KPfk_PEP = parameters[6]
		LPfk = parameters[7]
		KPfk_ADP_c = parameters[8]
		KPfk_AMP_a = parameters[9]
		nPfk = parameters[10]


		A=1+PEP/KPfk_PEP+ADP/KPfk_ADP_b+AMP/KPfk_AMP_b
		B=1+ADP/KPfk_ADP_a+AMP/KPfk_AMP_a

		rate=Pfk*kPfk_cat*ATP*F6P/((ATP+KPfk_ATP_s*(1+ADP/KPfk_ADP_c))*(F6P+KPfk_F6P_s*A/B)*(1+LPfk/pow(1+F6P*B/(KPfk_F6P_s*A),nPfk)))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Fbp_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Fbp"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KFbp_PEP','nFbp','KFbp_FBP','LFbp','kFbp_cat']
		self.reaction_substrate_names = ['PEP','FBP','Fbp']
		self.substrates = ['FBP']
		self.products = ['F6P']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		PEP = substrates[0]
		FBP = substrates[1]
		Fbp = substrates[2]


		KFbp_PEP = parameters[0]
		nFbp = parameters[1]
		KFbp_FBP = parameters[2]
		LFbp = parameters[3]
		kFbp_cat = parameters[4]



		rate=Fbp*kFbp_cat*FBP/KFbp_FBP*pow(1+FBP/KFbp_FBP,nFbp-1)/(pow(1+FBP/KFbp_FBP,nFbp)+LFbp/pow(1+PEP/KFbp_PEP,nFbp))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Fba_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Fba"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kFba_cat','VFba_blf','KFba_DHAP','KFba_eq','KFba_GAP','KFba_GAP_inh','KFba_FBP']
		self.reaction_substrate_names = ['GAP','Fba','FBP']
		self.substrates = ['FBP']
		self.products = ['GAP','DAP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		GAP = substrates[0]
		Fba = substrates[1]
		FBP = substrates[2]


		kFba_cat = parameters[0]
		VFba_blf = parameters[1]
		KFba_DHAP = parameters[2]
		KFba_eq = parameters[3]
		KFba_GAP = parameters[4]
		KFba_GAP_inh = parameters[5]
		KFba_FBP = parameters[6]



		rate=Fba*kFba_cat*(FBP-pow(GAP,2)/KFba_eq)/(KFba_FBP+FBP+KFba_GAP*GAP/(KFba_eq*VFba_blf)+KFba_DHAP*GAP/(KFba_eq*VFba_blf)+FBP*GAP/KFba_GAP_inh+pow(GAP,2)/(KFba_eq*VFba_blf))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Eno_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Eno"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KEno_Pga2','KEno_Pep','KEno_eq','vEno_max']
		self.reaction_substrate_names = ['PGA2','PEP']
		self.substrates = ['PGA2']
		self.products = ['PEP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		PGA2 = substrates[0]
		PEP = substrates[1]


		KEno_Pga2 = parameters[0]
		KEno_Pep = parameters[1]
		KEno_eq = parameters[2]
		vEno_max = parameters[3]



		rate=((vEno_max*(PGA2-PEP/KEno_eq))/KEno_Pga2)/(1+(PGA2/KEno_Pga2)+(PEP/KEno_Pep))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Tpi_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Tpi"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['K_TPI_Keq','K_TPI_KmGAP','K_TPI_KmDAP','v_TPI_Vmax']
		self.reaction_substrate_names = ['GAP','DAP']
		self.substrates = ['DAP']
		self.products = ['GAP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		GAP = substrates[0]
		DAP = substrates[1]


		K_TPI_Keq = parameters[0]
		K_TPI_KmGAP = parameters[1]
		K_TPI_KmDAP = parameters[2]
		v_TPI_Vmax = parameters[3]



		rate=v_TPI_Vmax*(DAP-GAP/K_TPI_Keq)/K_TPI_KmDAP/(1+DAP/K_TPI_KmDAP+GAP/K_TPI_KmGAP)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Pgk_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Pgk"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['K_PGK_KmPGA3','v_PGK_Vmax','K_PGK_Keq','K_PGK_KmBPG','K_PGK_KmADPMg','K_PGK_KmATPMg']
		self.reaction_substrate_names = ['BPG','ADP','PGA3','ATP']
		self.substrates = ['BPG','ADP']
		self.products = ['PGA3','ATP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		BPG = substrates[0]
		ADP = substrates[1]
		PGA3 = substrates[2]
		ATP = substrates[3]


		K_PGK_KmPGA3 = parameters[0]
		v_PGK_Vmax = parameters[1]
		K_PGK_Keq = parameters[2]
		K_PGK_KmBPG = parameters[3]
		K_PGK_KmADPMg = parameters[4]
		K_PGK_KmATPMg = parameters[5]



		rate=v_PGK_Vmax*(ADP*BPG-ATP*PGA3/K_PGK_Keq)/(K_PGK_KmADPMg*K_PGK_KmBPG)/(1+ADP/K_PGK_KmADPMg+BPG/K_PGK_KmBPG+ADP/K_PGK_KmADPMg*BPG/K_PGK_KmBPG+ATP/K_PGK_KmATPMg+PGA3/K_PGK_KmPGA3+ATP/K_PGK_KmATPMg*PGA3/K_PGK_KmPGA3)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Gpm_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Gpm"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['K_GPM_KmPGA2','K_GPM_KmPGA3','v_GPM_Vmax','K_GPM_Keq']
		self.reaction_substrate_names = ['PGA2','PGA3']
		self.substrates = ['PGA3']
		self.products = ['PGA2']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		PGA2 = substrates[0]
		PGA3 = substrates[1]


		K_GPM_KmPGA2 = parameters[0]
		K_GPM_KmPGA3 = parameters[1]
		v_GPM_Vmax = parameters[2]
		K_GPM_Keq = parameters[3]



		rate=v_GPM_Vmax*(PGA3-PGA2/K_GPM_Keq)/K_GPM_KmPGA3/(1+PGA3/K_GPM_KmPGA3+PGA2/K_GPM_KmPGA2)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Gdh_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Gdh"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['K_GDH_Keq','K_GDH_KmNAD','v_GDH_Vmax','K_GDH_KmP','K_GDH_KmNADH','K_GDH_KmGAP','K_GDH_KmBPG']
		self.reaction_substrate_names = ['NADH','BPG','NAD','P','GAP']
		self.substrates = ['GAP']
		self.products = ['BPG']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		NADH = substrates[0]
		BPG = substrates[1]
		NAD = substrates[2]
		P = substrates[3]
		GAP = substrates[4]


		K_GDH_Keq = parameters[0]
		K_GDH_KmNAD = parameters[1]
		v_GDH_Vmax = parameters[2]
		K_GDH_KmP = parameters[3]
		K_GDH_KmNADH = parameters[4]
		K_GDH_KmGAP = parameters[5]
		K_GDH_KmBPG = parameters[6]



		rate=v_GDH_Vmax*(P*GAP*NAD-BPG*NADH/K_GDH_Keq)/(K_GDH_KmP*K_GDH_KmGAP*K_GDH_KmNAD)/((1+P/K_GDH_KmP)*(1+GAP/K_GDH_KmGAP)*(1+NAD/K_GDH_KmNAD)+(1+BPG/K_GDH_KmBPG)*(1+NADH/K_GDH_KmNADH)-1)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Pyk_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Pyk"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kPyk_cat','KPyk_ADP','LPyk','KPyk_PEP','nPyk','KPyk_AMP','KPyk_ATP','KPyk_FBP']
		self.reaction_substrate_names = ['PEP','FBP','ATP','Pyk','ADP','AMP']
		self.substrates = ['PEP','ADP']
		self.products = ['PYR','ATP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		PEP = substrates[0]
		FBP = substrates[1]
		ATP = substrates[2]
		Pyk = substrates[3]
		ADP = substrates[4]
		AMP = substrates[5]


		kPyk_cat = parameters[0]
		KPyk_ADP = parameters[1]
		LPyk = parameters[2]
		KPyk_PEP = parameters[3]
		nPyk = parameters[4]
		KPyk_AMP = parameters[5]
		KPyk_ATP = parameters[6]
		KPyk_FBP = parameters[7]



		rate=Pyk*kPyk_cat*PEP*pow(PEP/KPyk_PEP+1,nPyk-1)*ADP/(KPyk_PEP*(LPyk*pow((1+ATP/KPyk_ATP)/(FBP/KPyk_FBP+AMP/KPyk_AMP+1),nPyk)+pow(PEP/KPyk_PEP+1,nPyk))*(ADP+KPyk_ADP))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Pps_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Pps"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['nPps','KPps_PYR','LPps','kPps_cat','KPps_PEP']
		self.reaction_substrate_names = ['Pps','PEP','PYR']
		self.substrates = ['PYR','ATP']
		self.products = ['PEP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		Pps = substrates[0]
		PEP = substrates[1]
		PYR = substrates[2]


		nPps = parameters[0]
		KPps_PYR = parameters[1]
		LPps = parameters[2]
		kPps_cat = parameters[3]
		KPps_PEP = parameters[4]



		rate=Pps*kPps_cat*PYR/KPps_PYR*pow(1+PYR/KPps_PYR,nPps-1)/(pow(1+PYR/KPps_PYR,nPps)+LPps*pow(1+PEP/KPps_PEP,nPps))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Pdh_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Pdh"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KPdh_NAD_m','KPdh_CoA_m','KPdh_i','KPdh_PYR_m','KPdh_AcCoA_m','kPdh_cat','KPdh_NADH_m']
		self.reaction_substrate_names = ['AcCoA','NADH','NAD','PYR','CoA','Pdh']
		self.substrates = ['PYR']
		self.products = ['AcCoA']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		AcCoA = substrates[0]
		NADH = substrates[1]
		NAD = substrates[2]
		PYR = substrates[3]
		CoA = substrates[4]
		Pdh = substrates[5]


		KPdh_NAD_m = parameters[0]
		KPdh_CoA_m = parameters[1]
		KPdh_i = parameters[2]
		KPdh_PYR_m = parameters[3]
		KPdh_AcCoA_m = parameters[4]
		kPdh_cat = parameters[5]
		KPdh_NADH_m = parameters[6]



		rate=Pdh*kPdh_cat*(1/(1+KPdh_i*NADH/NAD))*(PYR/KPdh_PYR_m)*(NAD/KPdh_NAD_m)*(CoA/KPdh_CoA_m)/((1+PYR/KPdh_PYR_m)*(1+NAD/KPdh_NAD_m+NADH/KPdh_NADH_m)*(1+CoA/KPdh_CoA_m+AcCoA/KPdh_AcCoA_m))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Pta_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Pta"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KPta_CoA_i','KPta_Pi_i','KPta_AcCoA_i','KPta_eq','vPta_max','KPta_AcP_i','KPta_Pi_m','KPta_AcP_m']
		self.reaction_substrate_names = ['AcCoA','Pi','CoA','AcP']
		self.substrates = ['AcCoA']
		self.products = ['AcP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		AcCoA = substrates[0]
		Pi = substrates[1]
		CoA = substrates[2]
		AcP = substrates[3]


		KPta_CoA_i = parameters[0]
		KPta_Pi_i = parameters[1]
		KPta_AcCoA_i = parameters[2]
		KPta_eq = parameters[3]
		vPta_max = parameters[4]
		KPta_AcP_i = parameters[5]
		KPta_Pi_m = parameters[6]
		KPta_AcP_m = parameters[7]



		rate=vPta_max*(1/(KPta_AcCoA_i*KPta_Pi_m))*(AcCoA*Pi-AcP*CoA/KPta_eq)/(1+AcCoA/KPta_AcCoA_i+Pi/KPta_Pi_i+AcP/KPta_AcP_i+CoA/KPta_CoA_i+AcCoA*Pi/(KPta_AcCoA_i*KPta_Pi_m)+AcP*CoA/(KPta_AcP_m*KPta_CoA_i))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Ack_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Ack"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['vAck_max','KAck_eq','KAck_ATP_m','KAck_ACE_m','KAck_AcP_m','KAck_ADP_m']
		self.reaction_substrate_names = ['ATP','ADP','ACEex','AcP']
		self.substrates = ['AcP','ADP']
		self.products = ['ATP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		ATP = substrates[0]
		ADP = substrates[1]
		ACEex = substrates[2]
		AcP = substrates[3]


		vAck_max = parameters[0]
		KAck_eq = parameters[1]
		KAck_ATP_m = parameters[2]
		KAck_ACE_m = parameters[3]
		KAck_AcP_m = parameters[4]
		KAck_ADP_m = parameters[5]



		rate=vAck_max*(1/(KAck_ADP_m*KAck_AcP_m))*(AcP*ADP-ACEex*ATP/KAck_eq)/((1+AcP/KAck_AcP_m+ACEex/KAck_ACE_m)*(1+ADP/KAck_ADP_m+ATP/KAck_ATP_m))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Ack_medium_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Ack_medium"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['vAck_max','rho','KAck_eq','KAck_ATP_m','KAck_ACE_m','KAck_AcP_m','KAck_ADP_m']
		self.reaction_substrate_names = ['ATP','X','AcP','ADP','ACEex']
		self.substrates = []
		self.products = ['ACEex']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		ATP = substrates[0]
		X = substrates[1]
		AcP = substrates[2]
		ADP = substrates[3]
		ACEex = substrates[4]


		vAck_max = parameters[0]
		rho = parameters[1]
		KAck_eq = parameters[2]
		KAck_ATP_m = parameters[3]
		KAck_ACE_m = parameters[4]
		KAck_AcP_m = parameters[5]
		KAck_ADP_m = parameters[6]



		vE_Ack=vAck_max*(1/(KAck_ADP_m*KAck_AcP_m))*(AcP*ADP-ACEex*ATP/KAck_eq)/((1+AcP/KAck_AcP_m+ACEex/KAck_ACE_m)*(1+ADP/KAck_ADP_m+ATP/KAck_ATP_m))
		rate=vE_Ack*X/rho
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Acs_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Acs"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KAcs_ACE','kAcs_cat']
		self.reaction_substrate_names = ['Acs','ACEex']
		self.substrates = ['ATP']
		self.products = ['AcCoA']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		Acs = substrates[0]
		ACEex = substrates[1]


		KAcs_ACE = parameters[0]
		kAcs_cat = parameters[1]



		rate=Acs*kAcs_cat*ACEex/(ACEex+KAcs_ACE)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Acs_medium_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Acs_medium"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KAcs_ACE','rho','kAcs_cat']
		self.reaction_substrate_names = ['Acs','X','ACEex']
		self.substrates = ['ACEex']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		Acs = substrates[0]
		X = substrates[1]
		ACEex = substrates[2]


		KAcs_ACE = parameters[0]
		rho = parameters[1]
		kAcs_cat = parameters[2]



		vE_Acs=Acs*kAcs_cat*ACEex/(ACEex+KAcs_ACE)
		rate=vE_Acs*X/rho
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Cs_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Cs"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KCs_AcCoA','KCs_OAA','kCs_cat','KCs_OAA_AcCoA','KCs_aKG']
		self.reaction_substrate_names = ['Cs','AcCoA','aKG','OAA']
		self.substrates = ['AcCoA','OAA']
		self.products = ['ICIT']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		Cs = substrates[0]
		AcCoA = substrates[1]
		aKG = substrates[2]
		OAA = substrates[3]


		KCs_AcCoA = parameters[0]
		KCs_OAA = parameters[1]
		kCs_cat = parameters[2]
		KCs_OAA_AcCoA = parameters[3]
		KCs_aKG = parameters[4]



		rate=Cs*kCs_cat*OAA*AcCoA/((1+aKG/KCs_aKG)*KCs_OAA_AcCoA*KCs_AcCoA+KCs_AcCoA*OAA+(1+aKG/KCs_aKG)*KCs_OAA*AcCoA+OAA*AcCoA)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Icdh_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Icdh"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KIcdh_ICIT','KIcdh_PEP','icdh_b','kIcdh_cat','LIcdh','nIcdh']
		self.reaction_substrate_names = ['Icdh','ICIT','PEP']
		self.substrates = ['ICIT']
		self.products = ['aKG']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		Icdh = substrates[0]
		ICIT = substrates[1]
		PEP = substrates[2]


		KIcdh_ICIT = parameters[0]
		KIcdh_PEP = parameters[1]
		icdh_b = parameters[2]
		kIcdh_cat = parameters[3]
		LIcdh = parameters[4]
		nIcdh = parameters[5]



		rate=(Icdh+icdh_b)*kIcdh_cat*ICIT/KIcdh_ICIT*pow(1+ICIT/KIcdh_ICIT,nIcdh-1)/(pow(1+ICIT/KIcdh_ICIT,nIcdh)+LIcdh*pow(1+PEP/KIcdh_PEP,nIcdh))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_akgdh_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_akgdh"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['Kakgdh_Z','Kakgdh_CoA_m','kakgdh_cat','Kakgdh_aKG_I','Kakgdh_SUC_I','Kakgdh_aKG_m','Kakgdh_NADH_I','Kakgdh_NAD_m']
		self.reaction_substrate_names = ['NADH','akgdh','NAD','CoA','aKG','SUC']
		self.substrates = ['aKG','ADP']
		self.products = ['SUC','ATP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		NADH = substrates[0]
		akgdh = substrates[1]
		NAD = substrates[2]
		CoA = substrates[3]
		aKG = substrates[4]
		SUC = substrates[5]


		Kakgdh_Z = parameters[0]
		Kakgdh_CoA_m = parameters[1]
		kakgdh_cat = parameters[2]
		Kakgdh_aKG_I = parameters[3]
		Kakgdh_SUC_I = parameters[4]
		Kakgdh_aKG_m = parameters[5]
		Kakgdh_NADH_I = parameters[6]
		Kakgdh_NAD_m = parameters[7]



		rate=akgdh*kakgdh_cat*aKG*CoA*NAD/(Kakgdh_NAD_m*aKG*CoA+Kakgdh_CoA_m*aKG*NAD+Kakgdh_aKG_m*CoA*NAD+aKG*CoA*NAD+Kakgdh_aKG_m*Kakgdh_Z*SUC*NADH/Kakgdh_SUC_I+Kakgdh_NAD_m*aKG*CoA*NADH/Kakgdh_NADH_I+Kakgdh_CoA_m*aKG*NAD*SUC/Kakgdh_SUC_I+Kakgdh_aKG_m*Kakgdh_Z*aKG*SUC*NADH/(Kakgdh_aKG_I*Kakgdh_SUC_I))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Sdh_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Sdh"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KSdh_SUC_m','KSdh_eq','kSdh1_cat','kSdh2_cat']
		self.reaction_substrate_names = ['FUM','Sdh','SUC']
		self.substrates = ['SUC']
		self.products = ['FUM']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		FUM = substrates[0]
		Sdh = substrates[1]
		SUC = substrates[2]


		KSdh_SUC_m = parameters[0]
		KSdh_eq = parameters[1]
		kSdh1_cat = parameters[2]
		kSdh2_cat = parameters[3]


		vSdh1_max=Sdh*kSdh1_cat
		vSdh2_max=Sdh*kSdh2_cat

		rate=vSdh1_max*vSdh2_max*(SUC-FUM/KSdh_eq)/(KSdh_SUC_m*vSdh2_max+vSdh2_max*SUC+vSdh1_max*FUM/KSdh_eq)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Fum_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Fum"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kFum2_cat','KFum_FUM_m','kFum1_cat','KFum_eq']
		self.reaction_substrate_names = ['Fum','MAL','FUM']
		self.substrates = ['FUM']
		self.products = ['MAL']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		Fum = substrates[0]
		MAL = substrates[1]
		FUM = substrates[2]


		kFum2_cat = parameters[0]
		KFum_FUM_m = parameters[1]
		kFum1_cat = parameters[2]
		KFum_eq = parameters[3]


		vFum1_max=Fum*kFum1_cat
		vFum2_max=Fum*kFum2_cat

		rate=vFum1_max*vFum2_max*(FUM-MAL/KFum_eq)/(KFum_FUM_m*vFum2_max+vFum2_max*FUM+vFum1_max*MAL/KFum_eq)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Mdh_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Mdh"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KMdh_MAL_m','kMdh1_cat','KMdh_NAD_II','KMdh_OAA_II','KMdh_OAA_I','KMdh_NAD_I','KMdh_NAD_m','kMdh2_cat','KMdh_NADH_I','KMdh_eq','KMdh_MAL_I','KMdh_NADH_m','KMdh_OAA_m']
		self.reaction_substrate_names = ['MAL','NADH','NAD','OAA','Mdh']
		self.substrates = ['MAL']
		self.products = ['OAA']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		MAL = substrates[0]
		NADH = substrates[1]
		NAD = substrates[2]
		OAA = substrates[3]
		Mdh = substrates[4]


		KMdh_MAL_m = parameters[0]
		kMdh1_cat = parameters[1]
		KMdh_NAD_II = parameters[2]
		KMdh_OAA_II = parameters[3]
		KMdh_OAA_I = parameters[4]
		KMdh_NAD_I = parameters[5]
		KMdh_NAD_m = parameters[6]
		kMdh2_cat = parameters[7]
		KMdh_NADH_I = parameters[8]
		KMdh_eq = parameters[9]
		KMdh_MAL_I = parameters[10]
		KMdh_NADH_m = parameters[11]
		KMdh_OAA_m = parameters[12]


		vMdh1_max=Mdh*kMdh1_cat
		vMdh2_max=Mdh*kMdh2_cat

		rate=vMdh1_max*vMdh2_max*(NAD*MAL-NADH*OAA/KMdh_eq)/(KMdh_NAD_I*KMdh_MAL_m*vMdh2_max+KMdh_MAL_m*vMdh2_max*NAD+KMdh_NAD_m*vMdh2_max*MAL+vMdh2_max*NAD*MAL+KMdh_OAA_m*vMdh1_max*NADH/KMdh_eq+KMdh_NADH_m*vMdh1_max*OAA/KMdh_eq+vMdh1_max*NADH*OAA/KMdh_eq+vMdh1_max*KMdh_OAA_m*NAD*NADH/(KMdh_eq*KMdh_NAD_I)+vMdh2_max*KMdh_NAD_m*MAL*OAA/KMdh_OAA_I+vMdh2_max*NAD*MAL*NADH/KMdh_NADH_I+vMdh1_max*MAL*NADH*OAA/(KMdh_eq*KMdh_MAL_I)+vMdh2_max*NAD*MAL*OAA/KMdh_OAA_II+vMdh1_max*NAD*NADH*OAA/(KMdh_NAD_II*KMdh_eq)+KMdh_NAD_I*vMdh2_max*NAD*MAL*NADH*OAA/(KMdh_NAD_II*KMdh_OAA_m*KMdh_NADH_I))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_MaeB_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_MaeB"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['nMaeB','KMaeB_MAL','LMaeB','KMaeB_cAMP','KMaeB_AcCoA','kMaeB_cat']
		self.reaction_substrate_names = ['cAMP','MAL','MaeB','AcCoA']
		self.substrates = ['MAL']
		self.products = ['PYR']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		MAL = substrates[1]
		MaeB = substrates[2]
		AcCoA = substrates[3]


		nMaeB = parameters[0]
		KMaeB_MAL = parameters[1]
		LMaeB = parameters[2]
		KMaeB_cAMP = parameters[3]
		KMaeB_AcCoA = parameters[4]
		kMaeB_cat = parameters[5]



		rate=MaeB*kMaeB_cat*MAL/KMaeB_MAL*pow(1+MAL/KMaeB_MAL,nMaeB-1)/(pow(1+MAL/KMaeB_MAL,nMaeB)+LMaeB*pow(1+AcCoA/KMaeB_AcCoA+cAMP/KMaeB_cAMP,nMaeB))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Pck_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Pck"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KPck_OAA','KPck_ADP_i','KPck_ATP_i','KPck_OAA_I','kPck_cat','KPck_PEP','KPck_PEP_i','KPck_ATP_I']
		self.reaction_substrate_names = ['PEP','OAA','ATP','Pck','ADP']
		self.substrates = ['OAA','ATP']
		self.products = ['PEP','ADP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		PEP = substrates[0]
		OAA = substrates[1]
		ATP = substrates[2]
		Pck = substrates[3]
		ADP = substrates[4]


		KPck_OAA = parameters[0]
		KPck_ADP_i = parameters[1]
		KPck_ATP_i = parameters[2]
		KPck_OAA_I = parameters[3]
		kPck_cat = parameters[4]
		KPck_PEP = parameters[5]
		KPck_PEP_i = parameters[6]
		KPck_ATP_I = parameters[7]



		rate=Pck*kPck_cat*OAA*ATP/ADP/(KPck_OAA*ATP/ADP+OAA*ATP/ADP+KPck_ATP_i*KPck_OAA/KPck_ADP_i+KPck_ATP_i*KPck_OAA/(KPck_PEP*KPck_ADP_i)*PEP+KPck_ATP_i*KPck_OAA/(KPck_PEP_i*KPck_ATP_I)*ATP/ADP*PEP+KPck_ATP_i*KPck_OAA/(KPck_ADP_i*KPck_OAA_I)*OAA)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Ppc_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Ppc"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KPpc_PEP','nPpc','kPpc_cat','LPpc','KPpc_FBP']
		self.reaction_substrate_names = ['PEP','Ppc','FBP']
		self.substrates = ['PEP']
		self.products = ['OAA']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		PEP = substrates[0]
		Ppc = substrates[1]
		FBP = substrates[2]


		KPpc_PEP = parameters[0]
		nPpc = parameters[1]
		kPpc_cat = parameters[2]
		LPpc = parameters[3]
		KPpc_FBP = parameters[4]



		rate=Ppc*kPpc_cat*PEP/KPpc_PEP*pow(1+PEP/KPpc_PEP,nPpc-1)/(pow(1+PEP/KPpc_PEP,nPpc)+LPpc/pow(1+FBP/KPpc_FBP,nPpc))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Icl_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Icl"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KIcl_3PG','KIcl_PEP','LIcl','nIcl','KIcl_ICIT','kIcl_cat','KIcl_aKG']
		self.reaction_substrate_names = ['PEP','Icl','ICIT','GAP','aKG']
		self.substrates = ['ICIT']
		self.products = ['SUC','GOX']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		PEP = substrates[0]
		Icl = substrates[1]
		ICIT = substrates[2]
		GAP = substrates[3]
		aKG = substrates[4]


		KIcl_3PG = parameters[0]
		KIcl_PEP = parameters[1]
		LIcl = parameters[2]
		nIcl = parameters[3]
		KIcl_ICIT = parameters[4]
		kIcl_cat = parameters[5]
		KIcl_aKG = parameters[6]



		rate=Icl*kIcl_cat*ICIT/KIcl_ICIT*pow(1+ICIT/KIcl_ICIT,nIcl-1)/(pow(1+ICIT/KIcl_ICIT,nIcl)+LIcl*pow(1+PEP/KIcl_PEP+GAP/KIcl_3PG+aKG/KIcl_aKG,nIcl))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Ms_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Ms"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KMs_GOX_AcCoA','kMs_cat','KMs_AcCoA','KMs_GOX']
		self.reaction_substrate_names = ['AcCoA','GOX','Ms']
		self.substrates = ['AcCoA','GOX']
		self.products = ['MAL']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		AcCoA = substrates[0]
		GOX = substrates[1]
		Ms = substrates[2]


		KMs_GOX_AcCoA = parameters[0]
		kMs_cat = parameters[1]
		KMs_AcCoA = parameters[2]
		KMs_GOX = parameters[3]



		rate=Ms*kMs_cat*GOX*AcCoA/(KMs_GOX_AcCoA*KMs_AcCoA+KMs_AcCoA*GOX+KMs_GOX*AcCoA+GOX*AcCoA)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_AceKki_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_AceKki"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['LAceK','kAceKki_cat','nAceK','KAceK_OAA','KAceK_Icdh']
		self.reaction_substrate_names = ['AceK','OAA','Icdh']
		self.substrates = ['ICDH','ATP']
		self.products = ['ICDHP','ADP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		AceK = substrates[0]
		OAA = substrates[1]
		Icdh = substrates[2]


		LAceK = parameters[0]
		kAceKki_cat = parameters[1]
		nAceK = parameters[2]
		KAceK_OAA = parameters[3]
		KAceK_Icdh = parameters[4]



		rate=AceK*kAceKki_cat*Icdh/KAceK_Icdh*pow(1+Icdh/KAceK_Icdh,nAceK-1)/(pow(1+Icdh/KAceK_Icdh,nAceK)+LAceK*pow(1+OAA/KAceK_OAA,nAceK))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_AceKph_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_AceKph"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KAceKph_OAA','KAceK_IcdhP','nAceKph','LAceKph','kAceKph_cat']
		self.reaction_substrate_names = ['AceK','OAA','IcdhP']
		self.substrates = ['ICDHP']
		self.products = ['ICDH']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		AceK = substrates[0]
		OAA = substrates[1]
		IcdhP = substrates[2]


		KAceKph_OAA = parameters[0]
		KAceK_IcdhP = parameters[1]
		nAceKph = parameters[2]
		LAceKph = parameters[3]
		kAceKph_cat = parameters[4]



		rate=AceK*kAceKph_cat*IcdhP/KAceK_IcdhP*pow(1+IcdhP/KAceK_IcdhP,nAceKph-1)/(pow(1+IcdhP/KAceK_IcdhP,nAceKph)+LAceKph/pow(1+OAA/KAceKph_OAA,nAceKph))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_G6pdh_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_G6pdh"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KG6pdh_G6P','vG6pdh_max','KG6pdh_NADPH_nadpinh','KG6pdh_NADPH_g6pinh','KG6pdh_NADP']
		self.reaction_substrate_names = ['NADP','NADPH','G6P']
		self.substrates = ['G6P']
		self.products = ['sixPGL']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		NADP = substrates[0]
		NADPH = substrates[1]
		G6P = substrates[2]


		KG6pdh_G6P = parameters[0]
		vG6pdh_max = parameters[1]
		KG6pdh_NADPH_nadpinh = parameters[2]
		KG6pdh_NADPH_g6pinh = parameters[3]
		KG6pdh_NADP = parameters[4]



		rate=vG6pdh_max*G6P*NADP/((G6P+KG6pdh_G6P)*(1+NADPH/KG6pdh_NADPH_g6pinh)*(KG6pdh_NADP*(1+NADPH/KG6pdh_NADPH_nadpinh)+NADP))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Pgl_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Pgl"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KPgl_6PG_m','vPgl_max','KPgl_h2','KPgl_6PGL_m','KPgl_h1','KPgl_eq']
		self.reaction_substrate_names = ['cAMP','sixPG','sixPGL','FBP','pH','EIIAP','PYR','ATP']
		self.substrates = ['sixPGL']
		self.products = ['sixPG']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		sixPG = substrates[1]
		sixPGL = substrates[2]
		FBP = substrates[3]
		pH = substrates[4]
		EIIAP = substrates[5]
		PYR = substrates[6]
		ATP = substrates[7]


		KPgl_6PG_m = parameters[0]
		vPgl_max = parameters[1]
		KPgl_h2 = parameters[2]
		KPgl_6PGL_m = parameters[3]
		KPgl_h1 = parameters[4]
		KPgl_eq = parameters[5]



		H=pow(10,-pH)*1e+3
		rate=vPgl_max*(sixPGL-sixPG/KPgl_eq)/((1+H/KPgl_h1+KPgl_h2/H)*(KPgl_6PGL_m+sixPGL+KPgl_6PGL_m/KPgl_6PG_m*sixPG))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Edd_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Edd"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['pH_Edd_m','KEdd_KDPG_m','KEdd_6PG_m','vEdd_max','KEdd_eq','pK_Edd']
		self.reaction_substrate_names = ['KDPG','sixPG','pH']
		self.substrates = ['sixPG']
		self.products = ['KDPG']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		KDPG = substrates[0]
		sixPG = substrates[1]
		pH = substrates[2]


		pH_Edd_m = parameters[0]
		KEdd_KDPG_m = parameters[1]
		KEdd_6PG_m = parameters[2]
		vEdd_max = parameters[3]
		KEdd_eq = parameters[4]
		pK_Edd = parameters[5]


		QEdd_pH=1+2*pow(10,pH_Edd_m-pK_Edd)/(1+pow(10,pH-pK_Edd)+pow(10,2*pH_Edd_m-pH-pK_Edd))

		rate=vEdd_max*QEdd_pH*(sixPG-KDPG/KEdd_eq)/(KEdd_6PG_m+sixPG+KEdd_6PG_m*KDPG/KEdd_KDPG_m)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Eda_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Eda"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['pK_Eda','KEda_GAP_m','vEda_max','KEda_KDPG_m','pH_Eda_m','KEda_PYR_m','KEda_eq']
		self.reaction_substrate_names = ['KDPG','GAP','PYR','pH']
		self.substrates = ['KDPG']
		self.products = ['GAP','PYR']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		KDPG = substrates[0]
		GAP = substrates[1]
		PYR = substrates[2]
		pH = substrates[3]


		pK_Eda = parameters[0]
		KEda_GAP_m = parameters[1]
		vEda_max = parameters[2]
		KEda_KDPG_m = parameters[3]
		pH_Eda_m = parameters[4]
		KEda_PYR_m = parameters[5]
		KEda_eq = parameters[6]


		QEda_pH=1+2*pow(10,pH_Eda_m-pK_Eda)/(1+pow(10,pH-pK_Eda)+pow(10,2*pH_Eda_m-pH-pK_Eda))

		rate=vEda_max*QEda_pH*(KDPG-GAP*PYR/KEda_eq)/(KEda_KDPG_m+KDPG+KEda_KDPG_m*(PYR/KEda_PYR_m+GAP/KEda_GAP_m+PYR*GAP/(KEda_PYR_m*KEda_GAP_m)))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_6Pgdh_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_6Pgdh"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['v6Pgdh_max','K6Pgdh_NADPH_inh','K6Pgdh_6PG','K6Pgdh_NADP','K6Pgdh_ATP_inh']
		self.reaction_substrate_names = ['ATP','NADPH','NADP','sixPG']
		self.substrates = ['sixPG']
		self.products = ['RU5P']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		ATP = substrates[0]
		NADPH = substrates[1]
		NADP = substrates[2]
		sixPG = substrates[3]


		v6Pgdh_max = parameters[0]
		K6Pgdh_NADPH_inh = parameters[1]
		K6Pgdh_6PG = parameters[2]
		K6Pgdh_NADP = parameters[3]
		K6Pgdh_ATP_inh = parameters[4]



		rate=v6Pgdh_max*sixPG*NADP/((sixPG+K6Pgdh_6PG)*(NADP+K6Pgdh_NADP*(1+NADPH/K6Pgdh_NADPH_inh)*(1+ATP/K6Pgdh_ATP_inh)))
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_R5pi_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_R5pi"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['vR5pi_max','KR5pi_eq']
		self.reaction_substrate_names = ['R5P','RU5P']
		self.substrates = ['RU5P']
		self.products = ['R5P']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		R5P = substrates[0]
		RU5P = substrates[1]


		vR5pi_max = parameters[0]
		KR5pi_eq = parameters[1]



		rate=vR5pi_max*(RU5P-R5P/KR5pi_eq)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Rpe_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Rpe"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KRpe_eq','vRpe_max']
		self.reaction_substrate_names = ['X5P','RU5P']
		self.substrates = ['RU5P']
		self.products = ['X5P']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		X5P = substrates[0]
		RU5P = substrates[1]


		KRpe_eq = parameters[0]
		vRpe_max = parameters[1]



		rate=vRpe_max*(RU5P-X5P/KRpe_eq)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Tkt1_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Tkt1"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KTkt1_eq','vTkt1_max']
		self.reaction_substrate_names = ['S7P','R5P','X5P','GAP']
		self.substrates = ['R5P','X5P']
		self.products = ['GAP','S7P']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		S7P = substrates[0]
		R5P = substrates[1]
		X5P = substrates[2]
		GAP = substrates[3]


		KTkt1_eq = parameters[0]
		vTkt1_max = parameters[1]



		rate=vTkt1_max*(R5P*X5P-S7P*GAP/KTkt1_eq)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Tkt2_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Tkt2"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KTkt2_eq','vTkt2_max']
		self.reaction_substrate_names = ['X5P','GAP','E4P','F6P']
		self.substrates = ['X5P','E4P']
		self.products = ['F6P','GAP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		X5P = substrates[0]
		GAP = substrates[1]
		E4P = substrates[2]
		F6P = substrates[3]


		KTkt2_eq = parameters[0]
		vTkt2_max = parameters[1]



		rate=vTkt2_max*(X5P*E4P-F6P*GAP/KTkt2_eq)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Tal_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Tal"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['vTal_max','KTal_eq']
		self.reaction_substrate_names = ['GAP','S7P','E4P','F6P']
		self.substrates = ['GAP','S7P']
		self.products = ['F6P','E4P']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		GAP = substrates[0]
		S7P = substrates[1]
		E4P = substrates[2]
		F6P = substrates[3]


		vTal_max = parameters[0]
		KTal_eq = parameters[1]



		rate=vTal_max*(GAP*S7P-E4P*F6P/KTal_eq)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_Cya_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_Cya"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KCya_EIIAP','vCya_max']
		self.reaction_substrate_names = ['EIIAP']
		self.substrates = ['ATP']
		self.products = ['cAMP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		EIIAP = substrates[0]


		KCya_EIIAP = parameters[0]
		vCya_max = parameters[1]



		rate=vCya_max*EIIAP/(EIIAP+KCya_EIIAP)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vE_cAMPdegr_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vE_cAMPdegr"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KcAMPdegr_cAMP','vcAMPdegr_max']
		self.reaction_substrate_names = ['cAMP']
		self.substrates = ['cAMP']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]


		KcAMPdegr_cAMP = parameters[0]
		vcAMPdegr_max = parameters[1]



		rate=vcAMPdegr_max*cAMP/(cAMP+KcAMPdegr_cAMP)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class OP_NADH_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "OP_NADH"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['Kakgdh_Z','KMdh_MAL_m','KPdh_AcCoA_m','KPdh_PYR_m','kPdh_cat','Kakgdh_NADH_I','KPdh_NAD_m','kMdh1_cat','K_GDH_Keq','Kakgdh_aKG_I','KMdh_OAA_II','KMdh_NAD_m','kMdh2_cat','Kakgdh_aKG_m','K_GDH_KmGAP','KMdh_eq','KPdh_NADH_m','KPdh_CoA_m','kakgdh_cat','v_GDH_Vmax','POratio','KMdh_NADH_I','Kakgdh_NAD_m','KMdh_MAL_I','KMdh_NADH_m','KMdh_NAD_II','K_GDH_KmNAD','KPdh_i','K_GDH_KmP','Kakgdh_SUC_I','K_GDH_KmNADH','KMdh_OAA_I','KMdh_NAD_I','Kakgdh_CoA_m','K_GDH_KmBPG','KMdh_OAA_m']
		self.reaction_substrate_names = ['AcCoA','NADH','MAL','CoA','akgdh','BPG','NAD','OAA','PYR','P','SUC','Mdh','GAP','aKG','Pdh']
		self.substrates = ['ADP']
		self.products = ['ATP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		AcCoA = substrates[0]
		NADH = substrates[1]
		MAL = substrates[2]
		CoA = substrates[3]
		akgdh = substrates[4]
		BPG = substrates[5]
		NAD = substrates[6]
		OAA = substrates[7]
		PYR = substrates[8]
		P = substrates[9]
		SUC = substrates[10]
		Mdh = substrates[11]
		GAP = substrates[12]
		aKG = substrates[13]
		Pdh = substrates[14]


		Kakgdh_Z = parameters[0]
		KMdh_MAL_m = parameters[1]
		KPdh_AcCoA_m = parameters[2]
		KPdh_PYR_m = parameters[3]
		kPdh_cat = parameters[4]
		Kakgdh_NADH_I = parameters[5]
		KPdh_NAD_m = parameters[6]
		kMdh1_cat = parameters[7]
		K_GDH_Keq = parameters[8]
		Kakgdh_aKG_I = parameters[9]
		KMdh_OAA_II = parameters[10]
		KMdh_NAD_m = parameters[11]
		kMdh2_cat = parameters[12]
		Kakgdh_aKG_m = parameters[13]
		K_GDH_KmGAP = parameters[14]
		KMdh_eq = parameters[15]
		KPdh_NADH_m = parameters[16]
		KPdh_CoA_m = parameters[17]
		kakgdh_cat = parameters[18]
		v_GDH_Vmax = parameters[19]
		POratio = parameters[20]
		KMdh_NADH_I = parameters[21]
		Kakgdh_NAD_m = parameters[22]
		KMdh_MAL_I = parameters[23]
		KMdh_NADH_m = parameters[24]
		KMdh_NAD_II = parameters[25]
		K_GDH_KmNAD = parameters[26]
		KPdh_i = parameters[27]
		K_GDH_KmP = parameters[28]
		Kakgdh_SUC_I = parameters[29]
		K_GDH_KmNADH = parameters[30]
		KMdh_OAA_I = parameters[31]
		KMdh_NAD_I = parameters[32]
		Kakgdh_CoA_m = parameters[33]
		K_GDH_KmBPG = parameters[34]
		KMdh_OAA_m = parameters[35]



		vE_Gdh=v_GDH_Vmax*(P*GAP*NAD-BPG*NADH/K_GDH_Keq)/(K_GDH_KmP*K_GDH_KmGAP*K_GDH_KmNAD)/((1+P/K_GDH_KmP)*(1+GAP/K_GDH_KmGAP)*(1+NAD/K_GDH_KmNAD)+(1+BPG/K_GDH_KmBPG)*(1+NADH/K_GDH_KmNADH)-1)
		vE_Pdh=Pdh*kPdh_cat*(1/(1+KPdh_i*NADH/NAD))*(PYR/KPdh_PYR_m)*(NAD/KPdh_NAD_m)*(CoA/KPdh_CoA_m)/((1+PYR/KPdh_PYR_m)*(1+NAD/KPdh_NAD_m+NADH/KPdh_NADH_m)*(1+CoA/KPdh_CoA_m+AcCoA/KPdh_AcCoA_m))
		vE_akgdh=akgdh*kakgdh_cat*aKG*CoA*NAD/(Kakgdh_NAD_m*aKG*CoA+Kakgdh_CoA_m*aKG*NAD+Kakgdh_aKG_m*CoA*NAD+aKG*CoA*NAD+Kakgdh_aKG_m*Kakgdh_Z*SUC*NADH/Kakgdh_SUC_I+Kakgdh_NAD_m*aKG*CoA*NADH/Kakgdh_NADH_I+Kakgdh_CoA_m*aKG*NAD*SUC/Kakgdh_SUC_I+Kakgdh_aKG_m*Kakgdh_Z*aKG*SUC*NADH/(Kakgdh_aKG_I*Kakgdh_SUC_I))
		vMdh1_max=Mdh*kMdh1_cat
		vMdh2_max=Mdh*kMdh2_cat
		vE_Mdh=vMdh1_max*vMdh2_max*(NAD*MAL-NADH*OAA/KMdh_eq)/(KMdh_NAD_I*KMdh_MAL_m*vMdh2_max+KMdh_MAL_m*vMdh2_max*NAD+KMdh_NAD_m*vMdh2_max*MAL+vMdh2_max*NAD*MAL+KMdh_OAA_m*vMdh1_max*NADH/KMdh_eq+KMdh_NADH_m*vMdh1_max*OAA/KMdh_eq+vMdh1_max*NADH*OAA/KMdh_eq+vMdh1_max*KMdh_OAA_m*NAD*NADH/(KMdh_eq*KMdh_NAD_I)+vMdh2_max*KMdh_NAD_m*MAL*OAA/KMdh_OAA_I+vMdh2_max*NAD*MAL*NADH/KMdh_NADH_I+vMdh1_max*MAL*NADH*OAA/(KMdh_eq*KMdh_MAL_I)+vMdh2_max*NAD*MAL*OAA/KMdh_OAA_II+vMdh1_max*NAD*NADH*OAA/(KMdh_NAD_II*KMdh_eq)+KMdh_NAD_I*vMdh2_max*NAD*MAL*NADH*OAA/(KMdh_NAD_II*KMdh_OAA_m*KMdh_NADH_I))
		rate=(vE_Gdh+vE_Pdh+vE_akgdh+vE_Mdh)*POratio
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class OP_FADH2_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "OP_FADH2"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kSdh2_cat','POratio_prime','KSdh_eq','KSdh_SUC_m','kSdh1_cat']
		self.reaction_substrate_names = ['FUM','Sdh','SUC']
		self.substrates = ['ADP']
		self.products = ['ATP']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		FUM = substrates[0]
		Sdh = substrates[1]
		SUC = substrates[2]


		kSdh2_cat = parameters[0]
		POratio_prime = parameters[1]
		KSdh_eq = parameters[2]
		KSdh_SUC_m = parameters[3]
		kSdh1_cat = parameters[4]



		vSdh1_max=Sdh*kSdh1_cat
		vSdh2_max=Sdh*kSdh2_cat
		vE_Sdh=vSdh1_max*vSdh2_max*(SUC-FUM/KSdh_eq)/(KSdh_SUC_m*vSdh2_max+vSdh2_max*SUC+vSdh1_max*FUM/KSdh_eq)
		rate=vE_Sdh*POratio_prime
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class mu_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "mu"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['ATP']
		self.substrates = []
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		ATP = substrates[0]


		kATP = parameters[0]



		rate=kATP*ATP
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vgrowth_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vgrowth"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','X','PYR']
		self.substrates = []
		self.products = ['X']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		X = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*X
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_glk_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_glk"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['vglk_Cra_unbound','nCraFBP','KCraFBP','kATP','vglk_Cra_bound','Cratotal','Kglk_Cra']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['Glk']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		vglk_Cra_unbound = parameters[0]
		nCraFBP = parameters[1]
		KCraFBP = parameters[2]
		kATP = parameters[3]
		vglk_Cra_bound = parameters[4]
		Cratotal = parameters[5]
		Kglk_Cra = parameters[6]



		mu=kATP*ATP
		mu=kATP*ATP
		CraFBP=Cratotal*pow(FBP,nCraFBP)/(pow(FBP,nCraFBP)+pow(KCraFBP,nCraFBP))
		Cra=Cratotal-CraFBP
		rate=mu*((1-Cra/(Cra+Kglk_Cra))*vglk_Cra_unbound+Cra/(Cra+Kglk_Cra)*vglk_Cra_bound)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_pfkA_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_pfkA"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KpfkA_Cra','nCraFBP','vpfkA_Cra_unbound','KCraFBP','kATP','vpfkA_Cra_bound','Cratotal']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['Pfk']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		KpfkA_Cra = parameters[0]
		nCraFBP = parameters[1]
		vpfkA_Cra_unbound = parameters[2]
		KCraFBP = parameters[3]
		kATP = parameters[4]
		vpfkA_Cra_bound = parameters[5]
		Cratotal = parameters[6]



		mu=kATP*ATP
		mu=kATP*ATP
		CraFBP=Cratotal*pow(FBP,nCraFBP)/(pow(FBP,nCraFBP)+pow(KCraFBP,nCraFBP))
		Cra=Cratotal-CraFBP
		rate=mu*((1-Cra/(Cra+KpfkA_Cra))*vpfkA_Cra_unbound+Cra/(Cra+KpfkA_Cra)*vpfkA_Cra_bound)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_fbp_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_fbp"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['vfbp_Cra_bound','nCraFBP','KCraFBP','kATP','vfbp_Cra_unbound','Cratotal','Kfbp_Cra']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['Fbp']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		vfbp_Cra_bound = parameters[0]
		nCraFBP = parameters[1]
		KCraFBP = parameters[2]
		kATP = parameters[3]
		vfbp_Cra_unbound = parameters[4]
		Cratotal = parameters[5]
		Kfbp_Cra = parameters[6]



		mu=kATP*ATP
		mu=kATP*ATP
		CraFBP=Cratotal*pow(FBP,nCraFBP)/(pow(FBP,nCraFBP)+pow(KCraFBP,nCraFBP))
		Cra=Cratotal-CraFBP
		rate=mu*((1-Cra/(Cra+Kfbp_Cra))*vfbp_Cra_unbound+Cra/(Cra+Kfbp_Cra)*vfbp_Cra_bound)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_fbaA_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_fbaA"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KfbaA_Cra','vfbaA_Crp_unbound','nCraFBP','KCraFBP','Crptotal','kATP','vfbaA_Crp_bound','vfbaA_Cra_bound','KfbaA_Crp','KCrpcAMP','nCrpcAMP','Cratotal','vfbaA_Cra_unbound']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['Fba']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		KfbaA_Cra = parameters[0]
		vfbaA_Crp_unbound = parameters[1]
		nCraFBP = parameters[2]
		KCraFBP = parameters[3]
		Crptotal = parameters[4]
		kATP = parameters[5]
		vfbaA_Crp_bound = parameters[6]
		vfbaA_Cra_bound = parameters[7]
		KfbaA_Crp = parameters[8]
		KCrpcAMP = parameters[9]
		nCrpcAMP = parameters[10]
		Cratotal = parameters[11]
		vfbaA_Cra_unbound = parameters[12]



		mu=kATP*ATP
		CrpcAMP=Crptotal*pow(cAMP,nCrpcAMP)/(pow(cAMP,nCrpcAMP)+pow(KCrpcAMP,nCrpcAMP))
		mu=kATP*ATP
		CraFBP=Cratotal*pow(FBP,nCraFBP)/(pow(FBP,nCraFBP)+pow(KCraFBP,nCraFBP))
		Cra=Cratotal-CraFBP
		rate=mu*((1-Cra/(Cra+KfbaA_Cra))*vfbaA_Cra_unbound+Cra/(Cra+KfbaA_Cra)*vfbaA_Cra_bound+(1-CrpcAMP/(CrpcAMP+KfbaA_Crp))*vfbaA_Crp_unbound+CrpcAMP/(CrpcAMP+KfbaA_Crp)*vfbaA_Crp_bound)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_gapA_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_gapA"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['vgapA_Cra_bound','KgapA_Cra','nCraFBP','vgapA_Cra_unbound','KCraFBP','vgapA_Crp_unbound','kATP','Crptotal','KgapA_Crp','KCrpcAMP','nCrpcAMP','Cratotal','vgapA_Crp_bound']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['GAPDH']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		vgapA_Cra_bound = parameters[0]
		KgapA_Cra = parameters[1]
		nCraFBP = parameters[2]
		vgapA_Cra_unbound = parameters[3]
		KCraFBP = parameters[4]
		vgapA_Crp_unbound = parameters[5]
		kATP = parameters[6]
		Crptotal = parameters[7]
		KgapA_Crp = parameters[8]
		KCrpcAMP = parameters[9]
		nCrpcAMP = parameters[10]
		Cratotal = parameters[11]
		vgapA_Crp_bound = parameters[12]



		mu=kATP*ATP
		CrpcAMP=Crptotal*pow(cAMP,nCrpcAMP)/(pow(cAMP,nCrpcAMP)+pow(KCrpcAMP,nCrpcAMP))
		mu=kATP*ATP
		CraFBP=Cratotal*pow(FBP,nCraFBP)/(pow(FBP,nCraFBP)+pow(KCraFBP,nCraFBP))
		Cra=Cratotal-CraFBP
		rate=mu*((1-Cra/(Cra+KgapA_Cra))*vgapA_Cra_unbound+Cra/(Cra+KgapA_Cra)*vgapA_Cra_bound+(1-CrpcAMP/(CrpcAMP+KgapA_Crp))*vgapA_Crp_unbound+CrpcAMP/(CrpcAMP+KgapA_Crp)*vgapA_Crp_bound)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_pykF_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_pykF"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['vpykF_Cra_bound','nCraFBP','KCraFBP','kATP','vpykF_Cra_unbound','Cratotal','KpykF_Cra']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['Pyk']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		vpykF_Cra_bound = parameters[0]
		nCraFBP = parameters[1]
		KCraFBP = parameters[2]
		kATP = parameters[3]
		vpykF_Cra_unbound = parameters[4]
		Cratotal = parameters[5]
		KpykF_Cra = parameters[6]



		mu=kATP*ATP
		mu=kATP*ATP
		CraFBP=Cratotal*pow(FBP,nCraFBP)/(pow(FBP,nCraFBP)+pow(KCraFBP,nCraFBP))
		Cra=Cratotal-CraFBP
		rate=mu*((1-Cra/(Cra+KpykF_Cra))*vpykF_Cra_unbound+Cra/(Cra+KpykF_Cra)*vpykF_Cra_bound)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_ppsA_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_ppsA"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['vppsA_Cra_bound','KppsA_Cra','nCraFBP','vppsA_Cra_unbound','KCraFBP','kATP','Cratotal']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['Pps']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		vppsA_Cra_bound = parameters[0]
		KppsA_Cra = parameters[1]
		nCraFBP = parameters[2]
		vppsA_Cra_unbound = parameters[3]
		KCraFBP = parameters[4]
		kATP = parameters[5]
		Cratotal = parameters[6]



		mu=kATP*ATP
		mu=kATP*ATP
		CraFBP=Cratotal*pow(FBP,nCraFBP)/(pow(FBP,nCraFBP)+pow(KCraFBP,nCraFBP))
		Cra=Cratotal-CraFBP
		rate=mu*((1-Cra/(Cra+KppsA_Cra))*vppsA_Cra_unbound+Cra/(Cra+KppsA_Cra)*vppsA_Cra_bound)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_lpd_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_lpd"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['vlpd_PdhR_bound','Klpd_PdhR','vlpd_PdhR_unbound','PdhRtotal','KPdhRPYR','kATP','nPdhRPYR']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['PDH']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		vlpd_PdhR_bound = parameters[0]
		Klpd_PdhR = parameters[1]
		vlpd_PdhR_unbound = parameters[2]
		PdhRtotal = parameters[3]
		KPdhRPYR = parameters[4]
		kATP = parameters[5]
		nPdhRPYR = parameters[6]



		mu=kATP*ATP
		PdhRPYR=PdhRtotal*pow(PYR,nPdhRPYR)/(pow(PYR,nPdhRPYR)+pow(KPdhRPYR,nPdhRPYR))
		PdhR=PdhRtotal-PdhRPYR
		mu=kATP*ATP
		rate=mu*((1-PdhR/(PdhR+Klpd_PdhR))*vlpd_PdhR_unbound+PdhR/(PdhR+Klpd_PdhR)*vlpd_PdhR_bound)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_acs_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_acs"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['nacs','vacs_Crp_bound','vacs_Crp_unbound','Crptotal','kATP','KCrpcAMP','Kacs_Crp','nCrpcAMP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['Acs']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		nacs = parameters[0]
		vacs_Crp_bound = parameters[1]
		vacs_Crp_unbound = parameters[2]
		Crptotal = parameters[3]
		kATP = parameters[4]
		KCrpcAMP = parameters[5]
		Kacs_Crp = parameters[6]
		nCrpcAMP = parameters[7]



		mu=kATP*ATP
		CrpcAMP=Crptotal*pow(cAMP,nCrpcAMP)/(pow(cAMP,nCrpcAMP)+pow(KCrpcAMP,nCrpcAMP))
		mu=kATP*ATP
		rate=mu*((1-pow(CrpcAMP,nacs)/(pow(CrpcAMP,nacs)+pow(Kacs_Crp,nacs)))*vacs_Crp_unbound+pow(CrpcAMP,nacs)/(pow(CrpcAMP,nacs)+pow(Kacs_Crp,nacs))*vacs_Crp_bound)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_gltA_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_gltA"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['ngltA','vgltA_Crp_bound','kATP','Crptotal','vgltA_Crp_unbound','KCrpcAMP','KgltA_Crp','nCrpcAMP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['CS']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		ngltA = parameters[0]
		vgltA_Crp_bound = parameters[1]
		kATP = parameters[2]
		Crptotal = parameters[3]
		vgltA_Crp_unbound = parameters[4]
		KCrpcAMP = parameters[5]
		KgltA_Crp = parameters[6]
		nCrpcAMP = parameters[7]



		mu=kATP*ATP
		CrpcAMP=Crptotal*pow(cAMP,nCrpcAMP)/(pow(cAMP,nCrpcAMP)+pow(KCrpcAMP,nCrpcAMP))
		mu=kATP*ATP
		rate=mu*((1-pow(CrpcAMP,ngltA)/(pow(CrpcAMP,ngltA)+pow(KgltA_Crp,ngltA)))*vgltA_Crp_unbound+pow(CrpcAMP,ngltA)/(pow(CrpcAMP,ngltA)+pow(KgltA_Crp,ngltA))*vgltA_Crp_bound)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_icdA_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_icdA"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['nCraFBP','vicdA_Cra_bound','vicdA_Cra_unbound','KCraFBP','KicdA_Cra','kATP','Cratotal']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['ICDH']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		nCraFBP = parameters[0]
		vicdA_Cra_bound = parameters[1]
		vicdA_Cra_unbound = parameters[2]
		KCraFBP = parameters[3]
		KicdA_Cra = parameters[4]
		kATP = parameters[5]
		Cratotal = parameters[6]



		mu=kATP*ATP
		mu=kATP*ATP
		CraFBP=Cratotal*pow(FBP,nCraFBP)/(pow(FBP,nCraFBP)+pow(KCraFBP,nCraFBP))
		Cra=Cratotal-CraFBP
		rate=mu*((1-Cra/(Cra+KicdA_Cra))*vicdA_Cra_unbound+Cra/(Cra+KicdA_Cra)*vicdA_Cra_bound)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_sucAB_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_sucAB"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['KsucAB_Crp','Crptotal','kATP','vsucAB_Crp_unbound','vsucAB_Crp_bound','KCrpcAMP','nCrpcAMP','nsucAB']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['aKGDH']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		KsucAB_Crp = parameters[0]
		Crptotal = parameters[1]
		kATP = parameters[2]
		vsucAB_Crp_unbound = parameters[3]
		vsucAB_Crp_bound = parameters[4]
		KCrpcAMP = parameters[5]
		nCrpcAMP = parameters[6]
		nsucAB = parameters[7]



		mu=kATP*ATP
		CrpcAMP=Crptotal*pow(cAMP,nCrpcAMP)/(pow(cAMP,nCrpcAMP)+pow(KCrpcAMP,nCrpcAMP))
		mu=kATP*ATP
		rate=mu*((1-pow(CrpcAMP,nsucAB)/(pow(CrpcAMP,nsucAB)+pow(KsucAB_Crp,nsucAB)))*vsucAB_Crp_unbound+pow(CrpcAMP,nsucAB)/(pow(CrpcAMP,nsucAB)+pow(KsucAB_Crp,nsucAB))*vsucAB_Crp_bound)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_sdhCDAB_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_sdhCDAB"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['vsdhCDAB_Crp_bound','KsdhCDAB_Crp','nsdhCDAB','Crptotal','kATP','vsdhCDAB_Crp_unbound','KCrpcAMP','nCrpcAMP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['SDH']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		vsdhCDAB_Crp_bound = parameters[0]
		KsdhCDAB_Crp = parameters[1]
		nsdhCDAB = parameters[2]
		Crptotal = parameters[3]
		kATP = parameters[4]
		vsdhCDAB_Crp_unbound = parameters[5]
		KCrpcAMP = parameters[6]
		nCrpcAMP = parameters[7]



		mu=kATP*ATP
		CrpcAMP=Crptotal*pow(cAMP,nCrpcAMP)/(pow(cAMP,nCrpcAMP)+pow(KCrpcAMP,nCrpcAMP))
		mu=kATP*ATP
		rate=mu*((1-pow(CrpcAMP,nsdhCDAB)/(pow(CrpcAMP,nsdhCDAB)+pow(KsdhCDAB_Crp,nsdhCDAB)))*vsdhCDAB_Crp_unbound+pow(CrpcAMP,nsdhCDAB)/(pow(CrpcAMP,nsdhCDAB)+pow(KsdhCDAB_Crp,nsdhCDAB))*vsdhCDAB_Crp_bound)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_fumABC_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_fumABC"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['vfumABC_Crp_bound','nfumABC','Crptotal','kATP','nCrpcAMP','KCrpcAMP','KfumABC_Crp','vfumABC_Crp_unbound']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['Fum']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		vfumABC_Crp_bound = parameters[0]
		nfumABC = parameters[1]
		Crptotal = parameters[2]
		kATP = parameters[3]
		nCrpcAMP = parameters[4]
		KCrpcAMP = parameters[5]
		KfumABC_Crp = parameters[6]
		vfumABC_Crp_unbound = parameters[7]



		mu=kATP*ATP
		CrpcAMP=Crptotal*pow(cAMP,nCrpcAMP)/(pow(cAMP,nCrpcAMP)+pow(KCrpcAMP,nCrpcAMP))
		mu=kATP*ATP
		rate=mu*((1-pow(CrpcAMP,nfumABC)/(pow(CrpcAMP,nfumABC)+pow(KfumABC_Crp,nfumABC)))*vfumABC_Crp_unbound+pow(CrpcAMP,nfumABC)/(pow(CrpcAMP,nfumABC)+pow(KfumABC_Crp,nfumABC))*vfumABC_Crp_bound)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_mdh_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_mdh"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['Kmdh_Crp','Crptotal','kATP','vmdh_Crp_bound','KCrpcAMP','nCrpcAMP','vmdh_Crp_unbound']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['MDH']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		Kmdh_Crp = parameters[0]
		Crptotal = parameters[1]
		kATP = parameters[2]
		vmdh_Crp_bound = parameters[3]
		KCrpcAMP = parameters[4]
		nCrpcAMP = parameters[5]
		vmdh_Crp_unbound = parameters[6]



		mu=kATP*ATP
		CrpcAMP=Crptotal*pow(cAMP,nCrpcAMP)/(pow(cAMP,nCrpcAMP)+pow(KCrpcAMP,nCrpcAMP))
		mu=kATP*ATP
		rate=mu*((1-CrpcAMP/(CrpcAMP+Kmdh_Crp))*vmdh_Crp_unbound+CrpcAMP/(CrpcAMP+Kmdh_Crp)*vmdh_Crp_bound)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_maeB_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_maeB"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['SS_MaeB','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['MaeB']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		SS_MaeB = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*SS_MaeB
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_pckA_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_pckA"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['nCraFBP','KCraFBP','kATP','vpckA_Cra_bound','vpckA_Cra_unbound','KpckA_Cra','Cratotal']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['Pck']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		nCraFBP = parameters[0]
		KCraFBP = parameters[1]
		kATP = parameters[2]
		vpckA_Cra_bound = parameters[3]
		vpckA_Cra_unbound = parameters[4]
		KpckA_Cra = parameters[5]
		Cratotal = parameters[6]



		mu=kATP*ATP
		mu=kATP*ATP
		CraFBP=Cratotal*pow(FBP,nCraFBP)/(pow(FBP,nCraFBP)+pow(KCraFBP,nCraFBP))
		Cra=Cratotal-CraFBP
		rate=mu*((1-Cra/(Cra+KpckA_Cra))*vpckA_Cra_unbound+Cra/(Cra+KpckA_Cra)*vpckA_Cra_bound)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_ppc_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_ppc"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['SS_Ppc','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = []
		self.products = ['Ppc']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		SS_Ppc = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*SS_Ppc
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_aceA_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_aceA"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['Kmu','KaceBAK_Crp','KCrpcAMP','nCrpcAMP','Cratotal','vaceBAK_Cra_unbound','n_mu','KCraFBP','KaceBAK_GOX','kaceBAK_cat_IclR','KaceBAK_DNA','LaceBAK','IclRtotal','KaceBAK_PYRprime','Crptotal','kATP','vaceBAK_Crp_unbound','nCraFBP','KaceBAK_Cra','vaceBAK_Cra_bound','KaceBAK_PYR','vaceBAK_Crp_bound']
		self.reaction_substrate_names = ['cAMP','aceBAK_DNA','PYR','pH','EIIAP','ATP','GOX','FBP']
		self.substrates = []
		self.products = ['Icl']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		aceBAK_DNA = substrates[1]
		PYR = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		GOX = substrates[6]
		FBP = substrates[7]


		Kmu = parameters[0]
		KaceBAK_Crp = parameters[1]
		KCrpcAMP = parameters[2]
		nCrpcAMP = parameters[3]
		Cratotal = parameters[4]
		vaceBAK_Cra_unbound = parameters[5]
		n_mu = parameters[6]
		KCraFBP = parameters[7]
		KaceBAK_GOX = parameters[8]
		kaceBAK_cat_IclR = parameters[9]
		KaceBAK_DNA = parameters[10]
		LaceBAK = parameters[11]
		IclRtotal = parameters[12]
		KaceBAK_PYRprime = parameters[13]
		Crptotal = parameters[14]
		kATP = parameters[15]
		vaceBAK_Crp_unbound = parameters[16]
		nCraFBP = parameters[17]
		KaceBAK_Cra = parameters[18]
		vaceBAK_Cra_bound = parameters[19]
		KaceBAK_PYR = parameters[20]
		vaceBAK_Crp_bound = parameters[21]



		mu=kATP*ATP
		mu=kATP*ATP
		CrpcAMP=Crptotal*pow(cAMP,nCrpcAMP)/(pow(cAMP,nCrpcAMP)+pow(KCrpcAMP,nCrpcAMP))
		CraFBP=Cratotal*pow(FBP,nCraFBP)/(pow(FBP,nCraFBP)+pow(KCraFBP,nCraFBP))
		Cra=Cratotal-CraFBP
		rate=mu*pow(Kmu,n_mu)/(pow(mu,n_mu)+pow(Kmu,n_mu))*((1-Cra/(Cra+KaceBAK_Cra))*vaceBAK_Cra_unbound+Cra/(Cra+KaceBAK_Cra)*vaceBAK_Cra_bound+(1-CrpcAMP/(CrpcAMP+KaceBAK_Crp))*vaceBAK_Crp_unbound+CrpcAMP/(CrpcAMP+KaceBAK_Crp)*vaceBAK_Crp_bound+(1-aceBAK_DNA/KaceBAK_DNA*(1+PYR/KaceBAK_PYRprime)/(1+1/LaceBAK*(GOX/KaceBAK_GOX)*(1+GOX/KaceBAK_GOX)+aceBAK_DNA/KaceBAK_DNA+PYR/KaceBAK_PYR+aceBAK_DNA/KaceBAK_DNA*PYR/KaceBAK_PYRprime))*kaceBAK_cat_IclR*IclRtotal)
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_aceB_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_aceB"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['Kmu','KaceBAK_Crp','KCrpcAMP','nCrpcAMP','Cratotal','vaceBAK_Cra_unbound','n_mu','KCraFBP','KaceBAK_GOX','kaceBAK_cat_IclR','KaceBAK_DNA','LaceBAK','IclRtotal','KaceBAK_PYRprime','Crptotal','kATP','vaceBAK_Crp_unbound','Factor_aceB','nCraFBP','KaceBAK_Cra','vaceBAK_Cra_bound','KaceBAK_PYR','vaceBAK_Crp_bound']
		self.reaction_substrate_names = ['cAMP','aceBAK_DNA','PYR','pH','EIIAP','ATP','GOX','FBP']
		self.substrates = []
		self.products = ['MS']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		aceBAK_DNA = substrates[1]
		PYR = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		GOX = substrates[6]
		FBP = substrates[7]


		Kmu = parameters[0]
		KaceBAK_Crp = parameters[1]
		KCrpcAMP = parameters[2]
		nCrpcAMP = parameters[3]
		Cratotal = parameters[4]
		vaceBAK_Cra_unbound = parameters[5]
		n_mu = parameters[6]
		KCraFBP = parameters[7]
		KaceBAK_GOX = parameters[8]
		kaceBAK_cat_IclR = parameters[9]
		KaceBAK_DNA = parameters[10]
		LaceBAK = parameters[11]
		IclRtotal = parameters[12]
		KaceBAK_PYRprime = parameters[13]
		Crptotal = parameters[14]
		kATP = parameters[15]
		vaceBAK_Crp_unbound = parameters[16]
		Factor_aceB = parameters[17]
		nCraFBP = parameters[18]
		KaceBAK_Cra = parameters[19]
		vaceBAK_Cra_bound = parameters[20]
		KaceBAK_PYR = parameters[21]
		vaceBAK_Crp_bound = parameters[22]



		mu=kATP*ATP
		CrpcAMP=Crptotal*pow(cAMP,nCrpcAMP)/(pow(cAMP,nCrpcAMP)+pow(KCrpcAMP,nCrpcAMP))
		CraFBP=Cratotal*pow(FBP,nCraFBP)/(pow(FBP,nCraFBP)+pow(KCraFBP,nCraFBP))
		Cra=Cratotal-CraFBP
		vG_aceA=mu*pow(Kmu,n_mu)/(pow(mu,n_mu)+pow(Kmu,n_mu))*((1-Cra/(Cra+KaceBAK_Cra))*vaceBAK_Cra_unbound+Cra/(Cra+KaceBAK_Cra)*vaceBAK_Cra_bound+(1-CrpcAMP/(CrpcAMP+KaceBAK_Crp))*vaceBAK_Crp_unbound+CrpcAMP/(CrpcAMP+KaceBAK_Crp)*vaceBAK_Crp_bound+(1-aceBAK_DNA/KaceBAK_DNA*(1+PYR/KaceBAK_PYRprime)/(1+1/LaceBAK*(GOX/KaceBAK_GOX)*(1+GOX/KaceBAK_GOX)+aceBAK_DNA/KaceBAK_DNA+PYR/KaceBAK_PYR+aceBAK_DNA/KaceBAK_DNA*PYR/KaceBAK_PYRprime))*kaceBAK_cat_IclR*IclRtotal)
		rate=Factor_aceB*vG_aceA
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vG_aceK_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vG_aceK"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['Kmu','KaceBAK_Crp','KCrpcAMP','nCrpcAMP','Cratotal','vaceBAK_Cra_unbound','n_mu','Factor_aceK','KCraFBP','KaceBAK_GOX','kaceBAK_cat_IclR','KaceBAK_DNA','LaceBAK','IclRtotal','KaceBAK_PYRprime','Crptotal','kATP','vaceBAK_Crp_unbound','nCraFBP','KaceBAK_Cra','vaceBAK_Cra_bound','KaceBAK_PYR','vaceBAK_Crp_bound']
		self.reaction_substrate_names = ['cAMP','aceBAK_DNA','PYR','pH','EIIAP','ATP','GOX','FBP']
		self.substrates = []
		self.products = ['AceK']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		aceBAK_DNA = substrates[1]
		PYR = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		GOX = substrates[6]
		FBP = substrates[7]


		Kmu = parameters[0]
		KaceBAK_Crp = parameters[1]
		KCrpcAMP = parameters[2]
		nCrpcAMP = parameters[3]
		Cratotal = parameters[4]
		vaceBAK_Cra_unbound = parameters[5]
		n_mu = parameters[6]
		Factor_aceK = parameters[7]
		KCraFBP = parameters[8]
		KaceBAK_GOX = parameters[9]
		kaceBAK_cat_IclR = parameters[10]
		KaceBAK_DNA = parameters[11]
		LaceBAK = parameters[12]
		IclRtotal = parameters[13]
		KaceBAK_PYRprime = parameters[14]
		Crptotal = parameters[15]
		kATP = parameters[16]
		vaceBAK_Crp_unbound = parameters[17]
		nCraFBP = parameters[18]
		KaceBAK_Cra = parameters[19]
		vaceBAK_Cra_bound = parameters[20]
		KaceBAK_PYR = parameters[21]
		vaceBAK_Crp_bound = parameters[22]



		mu=kATP*ATP
		CrpcAMP=Crptotal*pow(cAMP,nCrpcAMP)/(pow(cAMP,nCrpcAMP)+pow(KCrpcAMP,nCrpcAMP))
		CraFBP=Cratotal*pow(FBP,nCraFBP)/(pow(FBP,nCraFBP)+pow(KCraFBP,nCraFBP))
		Cra=Cratotal-CraFBP
		vG_aceA=mu*pow(Kmu,n_mu)/(pow(mu,n_mu)+pow(Kmu,n_mu))*((1-Cra/(Cra+KaceBAK_Cra))*vaceBAK_Cra_unbound+Cra/(Cra+KaceBAK_Cra)*vaceBAK_Cra_bound+(1-CrpcAMP/(CrpcAMP+KaceBAK_Crp))*vaceBAK_Crp_unbound+CrpcAMP/(CrpcAMP+KaceBAK_Crp)*vaceBAK_Crp_bound+(1-aceBAK_DNA/KaceBAK_DNA*(1+PYR/KaceBAK_PYRprime)/(1+1/LaceBAK*(GOX/KaceBAK_GOX)*(1+GOX/KaceBAK_GOX)+aceBAK_DNA/KaceBAK_DNA+PYR/KaceBAK_PYR+aceBAK_DNA/KaceBAK_DNA*PYR/KaceBAK_PYRprime))*kaceBAK_cat_IclR*IclRtotal)
		rate=Factor_aceK*vG_aceA
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_X_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_X"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['D']
		self.reaction_substrate_names = ['X']
		self.substrates = ['X']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		X = substrates[0]


		D = parameters[0]



		rate=D*X
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_GLCfeed_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_GLCfeed"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['D']
		self.reaction_substrate_names = ['GLCfeed']
		self.substrates = []
		self.products = ['GLCex']


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		GLCfeed = substrates[0]


		D = parameters[0]



		rate=D*GLCfeed
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_GLCex_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_GLCex"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['D']
		self.reaction_substrate_names = ['GLCex']
		self.substrates = ['GLCex']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		GLCex = substrates[0]


		D = parameters[0]



		rate=D*GLCex
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_GLC_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_GLC"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','GLC','ATP','PYR']
		self.substrates = ['GLC']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		GLC = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*GLC
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_G6P_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_G6P"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','G6P','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['G6P']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		G6P = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*G6P
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_F6P_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_F6P"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','F6P','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['F6P']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		F6P = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*F6P
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_FBP_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_FBP"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['FBP']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*FBP
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_GAP_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_GAP"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR','GAP']
		self.substrates = ['GAP']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]
		GAP = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*GAP
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_PEP_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_PEP"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','PEP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['PEP']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		PEP = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*PEP
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_PYR_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_PYR"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','PYR','pH','EIIAP','ATP','FBP']
		self.substrates = ['PYR']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		PYR = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		FBP = substrates[5]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*PYR
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_AcCoA_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_AcCoA"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','AcCoA','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['AcCoA']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		AcCoA = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*AcCoA
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_AcP_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_AcP"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR','AcP']
		self.substrates = ['AcP']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]
		AcP = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*AcP
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_ACEex_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_ACEex"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['D']
		self.reaction_substrate_names = ['ACEex']
		self.substrates = ['ACEex']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		ACEex = substrates[0]


		D = parameters[0]



		rate=D*ACEex
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_ICIT_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_ICIT"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR','ICIT']
		self.substrates = ['ICIT']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]
		ICIT = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*ICIT
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_aKG_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_aKG"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR','aKG']
		self.substrates = ['aKG']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]
		aKG = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*aKG
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_SUC_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_SUC"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR','SUC']
		self.substrates = ['SUC']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]
		SUC = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*SUC
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_FUM_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_FUM"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','FUM','PYR']
		self.substrates = ['FUM']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		FUM = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*FUM
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_MAL_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_MAL"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','MAL','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['MAL']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		MAL = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*MAL
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_OAA_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_OAA"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','OAA','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['OAA']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		OAA = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*OAA
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_GOX_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_GOX"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','GOX','PYR']
		self.substrates = ['GOX']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		GOX = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*GOX
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_6PGL_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_6PGL"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','sixPGL','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['sixPGL']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		sixPGL = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*sixPGL
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_6PG_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_6PG"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','sixPG','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['sixPG']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		sixPG = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*sixPG
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_KDPG_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_KDPG"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','KDPG','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['KDPG']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		KDPG = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*KDPG
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_RU5P_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_RU5P"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','RU5P','pH','EIIAP','ATP','FBP','PYR']
		self.substrates = ['RU5P']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		RU5P = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		FBP = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*RU5P
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_R5P_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_R5P"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','R5P','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['R5P']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		R5P = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*R5P
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_X5P_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_X5P"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR','X5P']
		self.substrates = ['X5P']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]
		X5P = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*X5P
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_S7P_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_S7P"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','S7P','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['S7P']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		S7P = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*S7P
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_E4P_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_E4P"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','E4P','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['E4P']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		E4P = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*E4P
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_cAMP_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_cAMP"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['cAMP']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]


		kATP = parameters[0]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=mu*cAMP
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Glk_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Glk"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','Glk','PYR']
		self.substrates = ['Glk']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		Glk = substrates[5]
		PYR = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Glk
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Pfk_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Pfk"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','Pfk','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['Pfk']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		Pfk = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Pfk
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Fbp_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Fbp"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR','Fbp']
		self.substrates = ['Fbp']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]
		Fbp = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Fbp
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Fba_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Fba"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR','Fba']
		self.substrates = ['Fba']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]
		Fba = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Fba
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Gapdh_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Gapdh"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR','Gapdh']
		self.substrates = ['GAPDH']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]
		Gapdh = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Gapdh
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Pyk_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Pyk"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR','Pyk']
		self.substrates = ['Pyk']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]
		Pyk = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Pyk
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Pps_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Pps"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','Pps','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['Pps']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		Pps = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Pps
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Pdh_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Pdh"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR','Pdh']
		self.substrates = ['PDH']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]
		Pdh = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Pdh
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Acs_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Acs"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','Acs','ATP','EIIAP','PYR']
		self.substrates = ['Acs']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		Acs = substrates[3]
		ATP = substrates[4]
		EIIAP = substrates[5]
		PYR = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Acs
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Cs_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Cs"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','Cs','ATP','PYR']
		self.substrates = ['CS']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		Cs = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Cs
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Icdh_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Icdh"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR','Icdh']
		self.substrates = ['ICDH']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]
		Icdh = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Icdh
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_IcdhP_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_IcdhP"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','IcdhP','pH','EIIAP','ATP','FBP','PYR']
		self.substrates = ['ICDHP']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		IcdhP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		FBP = substrates[5]
		PYR = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*IcdhP
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_akgdh_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_akgdh"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','akgdh','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['aKGDH']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		akgdh = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*akgdh
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Sdh_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Sdh"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','Sdh','PYR']
		self.substrates = ['SDH']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		Sdh = substrates[5]
		PYR = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Sdh
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Fum_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Fum"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','Fum','ATP','EIIAP','PYR']
		self.substrates = ['Fum']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		Fum = substrates[3]
		ATP = substrates[4]
		EIIAP = substrates[5]
		PYR = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Fum
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Mdh_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Mdh"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR','Mdh']
		self.substrates = ['MDH']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]
		Mdh = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Mdh
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_MaeB_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_MaeB"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','MaeB','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['MaeB']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		MaeB = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*MaeB
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Pck_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Pck"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR','Pck']
		self.substrates = ['Pck']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]
		Pck = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Pck
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Ppc_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Ppc"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','Ppc','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['Ppc']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		Ppc = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Ppc
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Icl_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Icl"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR','Icl']
		self.substrates = ['Icl']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]
		Icl = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Icl
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_Ms_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_Ms"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['cAMP','FBP','pH','EIIAP','ATP','PYR','Ms']
		self.substrates = ['MS']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		cAMP = substrates[0]
		FBP = substrates[1]
		pH = substrates[2]
		EIIAP = substrates[3]
		ATP = substrates[4]
		PYR = substrates[5]
		Ms = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*Ms
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vD_AceK_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vD_AceK"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kdegr','kATP']
		self.reaction_substrate_names = ['AceK','cAMP','FBP','pH','EIIAP','ATP','PYR']
		self.substrates = ['AceK']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		AceK = substrates[0]
		cAMP = substrates[1]
		FBP = substrates[2]
		pH = substrates[3]
		EIIAP = substrates[4]
		ATP = substrates[5]
		PYR = substrates[6]


		kdegr = parameters[0]
		kATP = parameters[1]



		mu=kATP*ATP
		mu=kATP*ATP
		rate=(mu+kdegr)*AceK
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vBM_G6P_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vBM_G6P"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kBM_GLC_G6P']
		self.reaction_substrate_names = ['G6P']
		self.substrates = ['G6P']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		G6P = substrates[0]


		kBM_GLC_G6P = parameters[0]



		rate=kBM_GLC_G6P*G6P
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vBM_F6P_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vBM_F6P"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kBM_GLC_F6P']
		self.reaction_substrate_names = ['F6P']
		self.substrates = ['F6P']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		F6P = substrates[0]


		kBM_GLC_F6P = parameters[0]



		rate=kBM_GLC_F6P*F6P
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vBM_GAP_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vBM_GAP"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kBM_GLC_GAP']
		self.reaction_substrate_names = ['GAP']
		self.substrates = ['GAP']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		GAP = substrates[0]


		kBM_GLC_GAP = parameters[0]



		rate=kBM_GLC_GAP*GAP
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vBM_PEP_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vBM_PEP"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kBM_GLC_PEP']
		self.reaction_substrate_names = ['PEP']
		self.substrates = ['PEP']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		PEP = substrates[0]


		kBM_GLC_PEP = parameters[0]



		rate=kBM_GLC_PEP*PEP
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vBM_PYR_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vBM_PYR"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kBM_GLC_PYR']
		self.reaction_substrate_names = ['PYR']
		self.substrates = ['PYR']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		PYR = substrates[0]


		kBM_GLC_PYR = parameters[0]



		rate=kBM_GLC_PYR*PYR
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vBM_AcCoA_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vBM_AcCoA"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kBM_GLC_AcCoA']
		self.reaction_substrate_names = ['AcCoA']
		self.substrates = ['AcCoA']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		AcCoA = substrates[0]


		kBM_GLC_AcCoA = parameters[0]



		rate=kBM_GLC_AcCoA*AcCoA
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vBM_aKG_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vBM_aKG"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kBM_GLC_aKG']
		self.reaction_substrate_names = ['aKG']
		self.substrates = ['aKG']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		aKG = substrates[0]


		kBM_GLC_aKG = parameters[0]



		rate=kBM_GLC_aKG*aKG
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vBM_SUC_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vBM_SUC"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kBM_GLC_SUC']
		self.reaction_substrate_names = ['SUC']
		self.substrates = ['SUC']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		SUC = substrates[0]


		kBM_GLC_SUC = parameters[0]



		rate=kBM_GLC_SUC*SUC
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vBM_FUM_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vBM_FUM"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kBM_GLC_FUM']
		self.reaction_substrate_names = ['FUM']
		self.substrates = ['FUM']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		FUM = substrates[0]


		kBM_GLC_FUM = parameters[0]



		rate=kBM_GLC_FUM*FUM
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vBM_OAA_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vBM_OAA"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kBM_GLC_OAA']
		self.reaction_substrate_names = ['OAA']
		self.substrates = ['OAA']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		OAA = substrates[0]


		kBM_GLC_OAA = parameters[0]



		rate=kBM_GLC_OAA*OAA
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vBM_R5P_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vBM_R5P"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kBM_GLC_R5P']
		self.reaction_substrate_names = ['R5P']
		self.substrates = ['R5P']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		R5P = substrates[0]


		kBM_GLC_R5P = parameters[0]



		rate=kBM_GLC_R5P*R5P
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

class vBM_E4P_Reaction(kinetics.Reaction):
	def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=""):
	
		super().__init__()
		if(metabolic_map_name==''):
			self.metabolic_map_name = "vBM_E4P"
		else:
			self.metabolic_map_name = metabolic_map_name

		self.last_time=0
		self.rates=[]
		self.parameter_names=['kBM_GLC_E4P']
		self.reaction_substrate_names = ['E4P']
		self.substrates = ['E4P']
		self.products = []


	def calculate_rate(self, substrates, parameters, rates_over_time, time):
		E4P = substrates[0]


		kBM_GLC_E4P = parameters[0]



		rate=kBM_GLC_E4P*E4P
		if(time!=self.last_time):
			self.last_time=time
			self.rates.append(rate)
			rates_over_time[self.metabolic_map_name]=self.rates
		return rate

