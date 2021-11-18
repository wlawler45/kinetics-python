import kinetics
from scipy.stats import reciprocal, uniform, norm

def parse_txt_file_parameters(filename):
    with open(filename, "r") as f:
        parameter_start=False
        parameter_dict={}
        inner_dict={}
        name=""
        for line in f:
            if("MODEL PARAMETERS" in line):
                parameter_start=True
                continue
            if("% WT" in line):
                parameter_start=False
            
            
            if(parameter_start):
                if("%%%" in line):
                    continue
                if("% Kinetic parameter for" in line):
                    if(name != ""):
                        parameter_dict[name]=inner_dict
                    inner_dict={}
                    #print(line)
                    splitter=line.split()
                    if(splitter[-1].isalnum() or "-" in splitter[-1]):
                        name=splitter[-1]
                        
                    
                    else:
                        for i in range(-1,-10,-1):
                            if(splitter[i].isalnum() or "-" in splitter[i]):
                                name=splitter[i]
                                
                                break
                    
                    print(name)
                else:
                    split2=line.split()
                    if(len(split2)>2 and not split2[-2].isalpha()):
                        inner_dict[split2[0]]=float(split2[-2])
                        #print(split2[0])
                        #print(split2[-2])
                    
                
            else:
                continue
            
    print(parameter_dict)
    
    
if __name__ == "__main__":  
    parse_txt_file_parameters("kuratadata.txt")