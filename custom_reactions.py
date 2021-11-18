import kinetics
import kinetics.Uncertainty as ua
import matplotlib.pyplot as plt
from scipy.stats import norm

class PGI_Reaction(kinetics.Reaction):

    def __init__(self,
        param1='', param2='', species1='', species2='',
        substrates=[], products=[]):

    # This is required to inherit from kinetics.Reaction
    super().__init__()

    # Set parameter and substrates names from the arguments passed in.  The order is important here.
    self.parameter_names=['vPgi_max', 'KPgi_eq','KPgi_G6P','KPgi_F6P','KPgi_F6P_6pginh','KPgi_G6P_6pginh']
    self.reaction_substrate_names = ['G6P', 'F6P','6PG']

    # Set the substrates and products from the arguments passed in.
    # Substrates are used up in the reaction, while produces are generated.
    self.substrates = ['G6P']
    self.products = ['F6P']

    def calculate_rate(self, substrates, parameters):

        # This function is used to calculate the rate at each time step in the model
        # It takes substrates and parameters as arguments, which are lists with the same order as we defined in __init__.

        # Substrates
        G6P = substrates[0]
        F6P = substrates[1]
        r6PG = substrates[2]

        # Parameters
        vPgi_max = parameters[0]
        KPgi_eq = parameters[1]
        KPgi_G6P = parameters[2]
        KPgi_F6P = parameters[3]
        Kpgi_F6P_6pginh =parameters[4]
        KPgi_G6P_6pginh=parameters[5]
        # This is where the rate equation goes.  An example is shown.
        rate = (vPgi_max*(G6P-(F6P/KPgi_eq)))/(Kpqi_g6p*(1+(F6P)/(KPgi_F6P*(1+r6PG/KPgi_F6P_6pginh))+r6PG/KPgi_G6P_6pginh)+G6P)

        return rate
        
class FBP_Reaction(kinetics.Reaction):

    def __init__(self,
        param1='', param2='', species1='', species2='',
        substrates=[], products=[]):

    # This is required to inherit from kinetics.Reaction
    super().__init__()

    # Set parameter and substrates names from the arguments passed in.  The order is important here.
    self.parameter_names=['kFbp_cat', 'KFbp_FBP','nFbp','LFbp','KFbp_PEP']
    self.reaction_substrate_names = ['Fbp', 'FBP','PEP']

    # Set the substrates and products from the arguments passed in.
    # Substrates are used up in the reaction, while produces are generated.
    self.substrates = ['FBP']
    self.products = ['F6P']

    def calculate_rate(self, substrates, parameters):

        # This function is used to calculate the rate at each time step in the model
        # It takes substrates and parameters as arguments, which are lists with the same order as we defined in __init__.

        # Substrates
        Fbp = substrates[0]
        FBP = substrates[1]
        PEP = substrates[2]

        # Parameters
        kFbp_cat = parameters[0]
        KFbp_FBP = parameters[1]
        nFbp = parameters[2]
        LFbp = parameters[3]
        KFbp_PEP =parameters[4]
        
        # This is where the rate equation goes.  An example is shown.
        rate = (Fbp*kFbp_cat*(FBP/KFbp_FBP)*(1+FBP/KFbp_FBP)**(nFbp-1))/((1+FBP/KFbp_FBP)**(nFbp) + LFbp/((1+FBP/KFbp_FBP)**(nFbp-1)))

        return rate
        
class FBP_enzyme_growth(kinetics.Reaction):

    def __init__(self,
        param1='', param2='', species1='', species2='',
        substrates=[], products=[]):

    # This is required to inherit from kinetics.Reaction
    super().__init__()

    # Set parameter and substrates names from the arguments passed in.  The order is important here.
    self.parameter_names=['mu','kexpr2','Kfbp_Cra', 'vfbp_Cra_unbound','vfbp_Cra_bound']
    self.reaction_substrate_names = ['Cra', 'Fbp']

    # Set the substrates and products from the arguments passed in.
    # Substrates are used up in the reaction, while produces are generated.
    self.substrates = ['Cra']
    self.products = ['Fbp']

    def calculate_rate(self, substrates, parameters):

        # This function is used to calculate the rate at each time step in the model
        # It takes substrates and parameters as arguments, which are lists with the same order as we defined in __init__.

        # Substrates
        Cra = substrates[0]
        Fbp = substrates[1]
        

        # Parameters
        mu = parameters[0]
        kexpr2 = parameters[1]
        Kfbp_Cra = parameters[2]
        vfbp_Cra_unbound = parameters[3]
        vfbp_Cra_bound =parameters[4]
        
        # This is where the rate equation goes.  An example is shown.
        rate = (mu*kexpr2*((1-Cra/(Cra+Kfbp_Cra))*vfbp_Cra_unbound+Cra/(Cra+Kfbp_Cra)*vfbp_Cra_bound))

        return rate
    
class PFK_reaction(kinetics.Reaction):

    def __init__(self,
        param1='', param2='', species1='', species2='',
        substrates=[], products=[]):

    # This is required to inherit from kinetics.Reaction
    super().__init__()

    # Set parameter and substrates names from the arguments passed in.  The order is important here.
    self.parameter_names=['kPfk_cat','KPfk_ATP_s','KPfk_ADP_c','KPfk_F6P_s', 'LPfk','nPfk','KPfk_PEP','KPfk_ADP_b','KPfk_AMP_b','KPfk_ADP_a','KPfk_AMP_a']
    self.reaction_substrate_names = ['Pfk', 'ATP','F6P','ADP','PEP','AMP']

    # Set the substrates and products from the arguments passed in.
    # Substrates are used up in the reaction, while produces are generated.
    self.substrates = ['F6P','ATP']
    self.products = ['FBP','ADP']

    def calculate_rate(self, substrates, parameters):

        # This function is used to calculate the rate at each time step in the model
        # It takes substrates and parameters as arguments, which are lists with the same order as we defined in __init__.

        # Substrates
        Pfk = substrates[0]
        ATP = substrates[1]
        F6P = substrates[2]
        ADP = substrates[3]
        PEP = substrates[4]
        AMP = substrates[5]
        
        

        # Parameters
        kPfk_cat = parameters[0]
        KPfk_ATP_s = parameters[1]
        KPfk_ADP_c = parameters[2]
        KPfk_F6P_s = parameters[3]
        LPfk =parameters[4]
        nPfk =parameters[5]
        KPfk_PEP=parameters[6]
        KPfk_ADP_b=parameters[7]
        KPfk_AMP_b=parameters[8]
        KPfk_ADP_a=parameters[9]
        KPfk_AMP_a=parameters[10]
        
        # This is where the rate equation goes.  An example is shown.
        A=1+PEP/KPfk_PEP +ADP/KPfk_ADP_b + AMP/KPfk_AMP_b
        B=1+ADP/KPfk_ADP_a +AMP/KPfk_AMP_a
        rate = (Pfk*kPfk_cat*ATP*F6P)/((ATP+KPfk_ATP_s*(1+ADP/KPfk_ADP_c))(F6P+KPfk_F6P_s*A/B)(1+LPfk/(1+F6P*B/(KPfk_F6P_s*A))))

        return rate
    
class PFK_enzyme_growth(kinetics.Reaction):

    def __init__(self,
        param1='', param2='', species1='', species2='',
        substrates=[], products=[]):

    # This is required to inherit from kinetics.Reaction
    super().__init__()

    # Set parameter and substrates names from the arguments passed in.  The order is important here.
    self.parameter_names=['muexpr2','KpfkA_Cra', 'vpfkA_Cra_unbound','vpfkA_Cra_bound']
    self.reaction_substrate_names = ['Cra', 'Pfk']

    # Set the substrates and products from the arguments passed in.
    # Substrates are used up in the reaction, while produces are generated.
    self.substrates = ['Cra']
    self.products = ['Pfk']

    def calculate_rate(self, substrates, parameters):

        # This function is used to calculate the rate at each time step in the model
        # It takes substrates and parameters as arguments, which are lists with the same order as we defined in __init__.

        # Substrates
        Cra = substrates[0]
        Pfk = substrates[1]
        

        # Parameters
        muexpr2 = parameters[0]
        KpfkA_Cra = parameters[1]
        vpfkA_Cra_unbound = parameters[2]
        vpfkA_Cra_bound =parameters[3]
        
        # This is where the rate equation goes.  An example is shown.
        rate = (muexpr2*((1-Cra/(Cra+KpfkA_Cra))*vpfkA_Cra_unbound+Cra/(Cra+KpfkA_Cra)*vpfkA_Cra_bound))

        return rate
        
        
class FBA_Reaction(kinetics.Reaction):

    def __init__(self,
        param1='', param2='', species1='', species2='',
        substrates=[], products=[]):

    # This is required to inherit from kinetics.Reaction
    super().__init__()

    # Set parameter and substrates names from the arguments passed in.  The order is important here.
    self.parameter_names=['kFba_cat', 'KFba_eq','KFba_FBP','KFba_GAP','VFba_blf','KFba_DHAP','KFba_GAP_inh']
    self.reaction_substrate_names = ['Fba', 'FBP','GAP']

    # Set the substrates and products from the arguments passed in.
    # Substrates are used up in the reaction, while produces are generated.
    self.substrates = ['FBP']
    self.products = ['GAP']

    def calculate_rate(self, substrates, parameters):

        # This function is used to calculate the rate at each time step in the model
        # It takes substrates and parameters as arguments, which are lists with the same order as we defined in __init__.

        # Substrates
        Fba = substrates[0]
        FBP = substrates[1]
        GAP = substrates[2]

        # Parameters
        kFba_cat = parameters[0]
        KFba_eq = parameters[1]
        KFba_FBP = parameters[2]
        KFba_GAP = parameters[3]
        VFba_blf = parameters[4]
        KFba_DHAP = parameters[5]
        KFba_GAP_inh = parameters[6]
        
        # This is where the rate equation goes.  An example is shown.
        rate = (Fba*kFba_cat*(FBP-(GAP**2 /KFba_eq)))/(KFba_FBP+FBP+(KFba_GAP*GAP)/(KFba_eq*VFba_blf)+(KFba_DHAP*GAP)/(KFba_eq*VFba_blf)+(FBP*GAP)/KFba_GAP_inh + (GAP**2)/(KFba_eq*VFba_blf))

        return rate
        
class PYK_Reaction(kinetics.Reaction):

    def __init__(self,
        param1='', param2='', species1='', species2='',
        substrates=[], products=[]):

    # This is required to inherit from kinetics.Reaction
    super().__init__()

    # Set parameter and substrates names from the arguments passed in.  The order is important here.
    self.parameter_names=['kPyk_cat', 'KPyk_PEP','nPyk','LPyk','KPyk_ATP','KPyk_FBP','KPyk_AMP','KPyk_ADP']
    self.reaction_substrate_names = ['Pyk', 'PEP','ADP','ATP','FBP','AMP']

    # Set the substrates and products from the arguments passed in.
    # Substrates are used up in the reaction, while produces are generated.
    self.substrates = ['PEP','ADP']
    self.products = ['PYR','ATP']

    def calculate_rate(self, substrates, parameters):

        # This function is used to calculate the rate at each time step in the model
        # It takes substrates and parameters as arguments, which are lists with the same order as we defined in __init__.

        # Substrates
        Pyk = substrates[0]
        PEP= substrates[1]
        ADP = substrates[2]
        ATP = substrates[3]
        FBP = substrates[4]
        AMP = substrates[5]

        # Parameters
        kPyk_cat = parameters[0]
        KPyk_PEP = parameters[1]
        nPyk = parameters[2]
        LPyk = parameters[3]
        KPyk_ATP = parameters[4]
        KPyk_FBP = parameters[5] 
        KPyk_AMP = parameters[6]
        KPyk_ADP = parameters[7]
        
        # This is where the rate equation goes.  An example is shown.
        rate = (Pyk*kPyk_cat*PEP*((PEP/KPyk_PEP) +1)**(nPyk-1) *ADP)/(KPyk_PEP(LPyk((1+ATP/KPyk_ATP)/(FBP/KPyk_FBP+AMP/KPyk_AMP +1))**nPyk +(PEP/KPyk_PEP +1)**nPyk)(ADP+KPyk_ADP))

        return rate