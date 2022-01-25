import escher
from escher import Builder
import cobra
from time import sleep

class map_generator():
    def __init__(self):
        self.builder = Builder()
        self.builder.map_name = 'iJO1366.Central metabolism'

        # some other things to try:
        self.builder.scroll_behavior = 'zoom'
        cobra_model = cobra.io.load_json_model('iJO1366.json')
        self.builder.model = cobra_model

        self.builder.hide_secondary_metabolites = True
        self.builder.hide_all_labels = False
        self.solution = self.builder.model.optimize()
        
        
    def draw_data(self,data,timestep):
        
        for i, flux_value in enumerate(self.solution.fluxes):
            self.solution.fluxes[i] = 0
        print("hello")
        print(data.keys())
        for key in data.keys():
            if(key in self.solution.fluxes.keys()):
                print(key)
                print(data[key][timestep])
                self.solution.fluxes[key] = data[key][timestep]
        # solution.fluxes['PGI'] = 10000
        # solution.fluxes['MDH'] = -10000
        # solution.fluxes['ENO'] = 5000

        # display fluxes on map
        self.builder.reaction_data = self.solution.fluxes

        self.builder.reaction_scale = [{'type': 'min', 'color': '#fa0000', 'size': 20}, 
                                  {'type': 'mean', 'color': '#000000', 'size': 20}, 
                                  {'type': 'max', 'color': '#15fa00', 'size': 20}]
        self.builder.reaction_styles = ['color', 'size', 'text']
        savestring="timestep_"+str(timestep)+".html"
        self.builder.save_html(savestring)
        
        
if __name__ == "__main__":
    genny=map_generator()