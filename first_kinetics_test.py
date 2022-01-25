import kinetics
from scipy.stats import reciprocal, uniform, norm
import re
import os

def reaction_class_creator_new(filename1, filename2,filename3,filename4,filename5):
    with open(filename1, "r") as f:
        model_states={}
        for line in f:
            if(line[0]=="%"):
                continue
            else:
                positive=[]
                negative=[]
                splitter=line.replace(";","").replace("_dot","").split()
                #print(splitter)
                for i in range(2,len(splitter)):
                    if('*' in splitter[i]):
                        dublesplit=splitter[i].split('*')
                        #print(dublesplit)
                        
                    if(splitter[i] != "-" and splitter[i] != "+"):
                        if('*' in splitter[i]):
                            dublesplit=splitter[i].split('*')
                            #print(dublesplit)
                            if(splitter[i-1]=="-"):
                                for i in range(int(dublesplit[0])):
                                    negative.append(dublesplit[1])
                            else:
                                for i in range(int(dublesplit[0])):
                                    positive.append(dublesplit[1])
                        else:
                            if(splitter[i-1]=="-"):
                                negative.append(splitter[i])
                            else:
                                positive.append(splitter[i])
                    model_states[splitter[0]]=[positive,negative]
                    #print(model_states)
                    
    with open(filename3,"r") as f:
        metabolite_list=[]
        for line in f:
            splitterlinger=line.split()
            metabolite_list.append(splitterlinger[0])
        #print(metabolite_list)
        
    with open(filename2, "r") as f:
        reaction_names=[]
        reviewed={}
        taggedones=[]
        
        for line in f:  
            output=[]
            splitting=line.replace(";","").split()
            if(line[0]=="%" or len(splitting)==0 or line[0]==" "):
                #print(line)
                continue
            
            reaction_names.append(splitting[0])
            preprocess=line.replace(";","").replace("power","").replace("(","").replace(")","").replace("\n","").replace("=","")
            finalized=re.split('[- + / , *]', preprocess)
            redefined=[]
            for value in finalized:
                if(not value.isnumeric()):
                    redefined.append(value)
                    
            
            results=set(redefined)
            results.remove('')
            #print(splitting[0])
            #print(results)
            results.remove(splitting[0])
            for item in results:
                output.append(item)
                
            reviewed[splitting[0]]=output
        taggedused=[]
        for key in reviewed.keys():
            for item in reviewed[key]:
                if(item in reaction_names):
                    taggedones.append(item)
                    taggedused.append(key)
        
        tagged=set(taggedones)
        tagged2=set(taggedused)
        #print(tagged)
        #print(tagged2)
        
    with open(filename5, "r") as f:
        model_variables={}
        model_equations={}
        model_variable_names=[]
        variable_metabs=[]
        
        for line in f:
            equations=[]
            mv_params=[]
            splitting=line.replace(";","").split()
            if(line[0]=="%" or len(splitting)==0 or line[0]==" "):
                #print(line)
                continue
            preprocess=line.replace(";","").replace("power","").replace("(","").replace(")","").replace("\n","").replace("=","")
            finalized=re.split('[- + / , *]', preprocess)
            redefined=[]
            for value in finalized:
                if(value in metabolite_list):
                    variable_metabs.append(value)
                if(not value.isnumeric() and value not in metabolite_list):
                    redefined.append(value)
                    
            for entry in redefined:
                if(entry in model_variable_names):
                    #print(redefined)
                    redefined+=model_variables[entry][0]
                    redefined.remove(entry)
                    equations.append(model_equations[entry])
                    #print(redefined)
            results=set(redefined)
            results.remove('')
            
            results.remove(splitting[0])
            mv_params=list(results)
            
            #print(mv_params)        
            if('1e' in mv_params):
                mv_params.remove('1e')
            #print(mv_params)
            equation=line.replace(";","").replace("power","pow").replace("\n","").replace(" ","")
            equations.append(equation)
            #print(equation)
            model_equations[splitting[0]]=equation
            model_variables[splitting[0]]=[mv_params,equations,variable_metabs]
            model_variable_names.append(splitting[0])
            
        #print(model_variables)
            
    
    returner={}        
    with open(filename2, "r") as f:
        with open(filename4, "w") as v:
            v.write("import kinetics\n")
            #v.write("from kinetics import *\n")
            #v.write("import kinetics.Uncertainty as ua\n")
            v.write("import matplotlib.pyplot as plt\n")
            v.write("from scipy.stats import norm\n\n")
            
            
            
            previous_equation1=None
            previous_equation2=None
            signal_int=0
            for line in f:
                equation_multiples=[]
                splitting=line.replace(";","").replace("power","pow").split()
                #print(splitting)
                if(line[0]=="%" or len(splitting)==0):
                    #print(line)
                    continue
                elif(line[0]==" "):
                    signal_int+=1
                    if(signal_int==1):
                        previous_equation1=line
                    else:
                        
                        previous_equation2=line
                    
                else:
                    v.write("class %s_Reaction(kinetics.Reaction):\n\t"%splitting[0])
                    #v.write("class %s_Reaction(Reaction):\n\t"%splitting[0])
                    # if(splitting[0] in tagged or splitting[0] in tagged2):
                        # v.write("def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],tempdata=''):\n\t")
                    # else:
                    #v.write("def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[]):\n\t")
                    v.write("def __init__(self,param1='',param2='',species1='',species2='',substrates=[],products=[],metabolic_map_name=\"\"):\n\t")
                    v.write("\n\t\tsuper().__init__()\n")
                    v.write("\t\tif(metabolic_map_name==''):\n")
                    v.write("\t\t\tself.metabolic_map_name = \"%s\"\n"%splitting[0])
                    v.write("\t\telse:\n")
                    v.write("\t\t\tself.metabolic_map_name = metabolic_map_name\n\n")
                    v.write("\t\tself.last_time=0\n")
                    v.write("\t\tself.rates=[]\n")
                    values=[]
                    # if(splitting[0] in tagged or splitting[0] in tagged2):
                        # v.write("\t\tself.tempdata=tempdata\n")
                    
                    if(signal_int==1):
                        preprocess1=previous_equation1.replace(";","").replace("power","").replace("(","").replace(")","").replace("\n","").replace("=","")
                        finalized1=re.split('[- + / , *]', preprocess1)
                        redefined1=[]
                        for value in finalized1:
                            if(not value.isnumeric()):
                                redefined1.append(value)
                                
                        values+=redefined1
                        
                        
                    elif(signal_int==2):
                        preprocess1=previous_equation1.replace(";","").replace("power","").replace("(","").replace(")","").replace("\n","").replace("=","")
                        finalized1=re.split('[- + / , *]', preprocess1)
                        redefined1=[]
                        for value in finalized1:
                            if(not value.isnumeric()):
                                redefined1.append(value)
                                
                        values+=redefined1
                        preprocess2=previous_equation2.replace(";","").replace("power","").replace("(","").replace(")","").replace("\n","").replace("=","")
                        finalized2=re.split('[- + / , *]', preprocess2)
                        redefined2=[]
                        for value in finalized2:
                            if(not value.isnumeric()):
                                redefined2.append(value)
                                
                        values+=redefined2
                        #print("hello")
                    
                    preprocess=line.replace(";","").replace("power","").replace("(","").replace(")","").replace("\n","").replace("=","")
                    finalized=re.split('[- + / , *]', preprocess)
                    redefined=[]
                    for value in finalized:
                        if(not value.isnumeric()):
                            redefined.append(value)
                            
                    values+=redefined
                    results=set(values)
                    results.remove('')
                    results.remove(splitting[0])
                    metabolites=[]
                    parameters=[]
                    for item in results:
                        if(item in metabolite_list):
                            metabolites.append(item)
                        else:
                            parameters.append(item)
                    
                    if(signal_int==1):
                        pre=previous_equation1.replace(";","").replace("power","pow").replace("\n","")
                        split1=pre.split()
                        
                        
                        if(split1[0] in parameters):
                            parameters.remove(split1[0])
                            #print("HIIII")
                            #print("%s found"%split1[0])
                        
                    if(signal_int==2):
                        pre=previous_equation1.replace(";","").replace("power","pow").replace("\n","")
                        pre2=previous_equation2.replace(";","").replace("power","pow").replace("\n","")
                        #print(pre)
                        #print(pre2)
                        
                        split1=pre.split()
                        split2=pre2.split()
                        #rint(split1[0])
                        #print(split2[0])
                        if(split1[0] in parameters):
                            parameters.remove(split1[0])
                            #print("%s found"%split1[0])
                        
                        if(split2[0] in parameters):
                            parameters.remove(split2[0])
                            #print(parameters)
                            #print("%s found"%split2[0])
                    #print(parameters)
                    v.write("\t\tself.parameter_names=[")
                    param_copy=parameters.copy()
                    param_temp=[]
                    subreact_names=[]
                    sub_vars=[]
                    for item in parameters:
                        if(item not in reaction_names):
                            param_temp.append(item)
                            #print(param_copy)
                        else:
                            subreact_names.append(item)
                        if(item in model_variable_names):
                            sub_vars.append(item)
                    #print(subreact_names)
                    #print(sub_vars)
                    parameters=param_temp
                    #print(subreact_names)
                    
                    if(len(subreact_names)!=0):
                        for sub in subreact_names:
                            
                            parameters+=returner[sub][1]
                            metabolites+=returner[sub][0]
                            #print(sub)
                            #print(returner[sub][1])
                    if(len(sub_vars)!=0):
                        for sub in sub_vars:
                            parameters+=model_variables[sub][0]
                            metabolites+=model_variables[sub][2]
                            #print(splitting[0])
                            #print(model_variables[sub][0])
                    param_setty=set(parameters)
                    parameter_new=[]
                    for item in param_setty:
                        if(item not in model_variable_names):
                            parameter_new.append(item)
                    parameters=parameter_new
                    
                    if(len(parameters)==0):
                        v.write("]\n")
                    else:
                        for i in range(len(parameters)-1):
                            v.write("'%s',"%parameters[i])
                        v.write("'%s']\n"%parameters[-1])
                    v.write("\t\tself.reaction_substrate_names = [")
                    metaset=set(metabolites)
                    metabolites=list(metaset)
                    
                    if(len(metabolites)==0):
                        v.write("]\n")
                    else:
                        for i in range(len(metabolites)-1):
                            v.write("'%s',"%metabolites[i])
                        v.write("'%s']\n"%metabolites[-1])
                    substrates=[]
                    products=[]
                    for key in model_states.keys():
                        if(splitting[0] in model_states[key][0]):
                            countprod=model_states[key][0].count(splitting[0])
                            for i in range(countprod):
                                products.append(key)
                        if(splitting[0] in model_states[key][1]):
                            countsub=model_states[key][1].count(splitting[0])
                            for i in range(countsub):
                                substrates.append(key)
                    
                    v.write("\t\tself.substrates = [")
                    if(len(substrates)==0):
                        v.write("]\n")
                    else:
                        for i in range(len(substrates)-1):
                            v.write("'%s',"%substrates[i])
                        v.write("'%s']\n"%substrates[-1])
                        
                        
                    v.write("\t\tself.products = [")
                    if(len(products)==0):
                        v.write("]\n")
                    else:
                        for i in range(len(products)-1):
                            v.write("'%s',"%products[i])
                        v.write("'%s']\n"%products[-1])
                        
                        
                    v.write("\n\n\tdef calculate_rate(self, substrates, parameters, rates_over_time, time):\n")
                    for i in range(len(metabolites)):
                        v.write("\t\t%s = substrates[%i]\n"%(metabolites[i],i))
                    v.write("\n\n")
                    for i in range(len(parameters)):
                        v.write("\t\t%s = parameters[%i]\n"%(parameters[i],i))
                    
                    # for item in param_copy:
                        # if(item in tagged):
                            # v.write("\t\t%s = self.tempdata['%s']\n"%(item,item))
                    
                    v.write("\n\n")
                    if(signal_int==1):
                        pre=previous_equation1.replace(";","").replace("power","pow").replace("\n","")
                        split1=pre.split()
                        prel=pre.replace(" ","")
                        v.write("\t\t%s\n"%prel)
                        equation_multiples.append(prel)
                        #print(split1[0])
                        if(split1[0] in parameters):
                            parameters.remove(split1[0])
                            #print("HIIII")
                            #print("%s found"%split1[0])
                        
                    if(signal_int==2):
                        pre=previous_equation1.replace(";","").replace("power","pow").replace("\n","")
                        pre2=previous_equation2.replace(";","").replace("power","pow").replace("\n","")
                        #print(pre)
                        #print(pre2)
                        prel=pre.replace(" ","")
                        prel2=pre2.replace(" ","")
                        v.write("\t\t%s\n"%prel)
                        v.write("\t\t%s\n"%prel2)
                        equation_multiples.append(prel)
                        equation_multiples.append(prel2)
                        split1=pre.split()
                        split2=pre2.split()
                        #rint(split1[0])
                        #print(split2[0])
                        if(split1[0] in parameters):
                            parameters.remove(split1[0])
                            #print("%s found"%split1[0])
                        if(split2[0] in parameters):
                            parameters.remove(split2[0])
                            #print(parameters)
                            #print("%s found"%split2[0])
                    
                    
                    equation=line.replace(";","").replace("power","pow").replace("\n","").replace(" ","").replace(splitting[0],"rate",1)
                    #print(equation)
                    if(len(subreact_names)!=0):
                        for sub in subreact_names:
                            if(len(returner[sub][4])!=0):
                                for val in returner[sub][4]:
                                    new_equals=model_variables[val][1]
                                    for equalize in new_equals:
                                        #print(equalize)
                                        v.write("\n\t\t%s"%equalize)
                            
                            for my_equation in returner[sub][2]:
                                redone_equation=my_equation.replace("rate",sub)
                            
                                #print(redone_equation)
                                v.write("\n\t\t%s"%redone_equation)
                            
                            #print(sub)
                            #print(returner[sub][1])
                    if(len(sub_vars)!=0):
                        for sub in sub_vars:
                            new_equals=model_variables[sub][1]
                            for equalize in new_equals:
                                #print(equalize)
                                v.write("\n\t\t%s"%equalize)
                            #print(splitting[0])
                            
                    v.write("\n\t\t%s\n"%equation)
                    equation_multiples.append(equation)
                    # if(splitting[0] in tagged):
                        # v.write("\n\t\tself.tempdata['%s']=rate\n"%(splitting[0]))
                    v.write("\t\tif(time!=self.last_time):\n")
                    v.write("\t\t\tself.last_time=time\n")
                    v.write("\t\t\tself.rates.append(rate)\n")
                    v.write("\t\t\trates_over_time[self.metabolic_map_name]=self.rates\n")
                    v.write("\t\treturn rate\n\n")
                    returner[splitting[0]]=[metabolites,parameters,equation_multiples,subreact_names,sub_vars]
                    #print(results)
                    
                    #print(results)
                    signal_int=0
    union=set.union(tagged2,tagged)                
    return returner,union
                
                


def parse_txt_file_parameters(filename):
    with open(filename, "r") as f:
        parameter_start=False
        parameter_dict={}
        inner_dict={}
        unidict={}
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
                    
                    #print(name)
                else:
                    split2=line.split()
                    if(len(split2)>2 and not split2[-2].isalpha()):
                        #print(line)
                        #TODO: Make optional inputs for units of kcat and vmax
                        if("cat" in split2[0]):
                            print(split2[0])
                            inner_dict[split2[0]]=float(split2[-2])/60.0
                            unidict[split2[0]]=float(split2[-2])/60.0
                        elif(split2[0][0]=='k' or split2[0][0]=='K'):
                            if("degr" in split2[0]):
                                print(split2[0])
                                inner_dict[split2[0]]=float(split2[-2])/60.0
                                unidict[split2[0]]=float(split2[-2])/60.0
                            elif("kBM" in split2[0]):
                                print(split2[0])
                                inner_dict[split2[0]]=float(split2[-2])/60.0
                                unidict[split2[0]]=float(split2[-2])/60.0
                            else:
                                inner_dict[split2[0]]=float(split2[-2])*1000.0
                                unidict[split2[0]]=float(split2[-2])*1000.0
                        elif("max" in split2[0]):
                            inner_dict[split2[0]]=float(split2[-2])/60.0
                            unidict[split2[0]]=float(split2[-2])/60.0
                        elif("bound" in split2[0]):
                            inner_dict[split2[0]]=(float(split2[-2])/60.0)*1000.0
                            unidict[split2[0]]=(float(split2[-2])/60.0)*1000.0
                        else:
                            inner_dict[split2[0]]=float(split2[-2])
                            unidict[split2[0]]=float(split2[-2])
                        
                        #print(split2[0])
                        #print(split2[-2])
                    
                
            else:
                continue
            
    return parameter_dict,unidict

def get_initial_conditions(filename):
    with open(filename, "r") as v:
        initials={}
        for line in v:
            splitter=line.split()
            initials[splitter[-1]]=float(splitter[2])*1000
        #print(initials)
        return initials
        
def get_mapping_names(filename):
    with open(filename, "r") as v:
        initials={}
        for line in v:
            splitter=line.split()
            initials[splitter[0]]=splitter[1]
        #print(initials)
        return initials


def write_running_code(filename,mapping_names,parameter_dict,returned,init,tagged):
    with open(filename, "w") as v:
        v.write("import kinetics\n")
        v.write("import matplotlib.pyplot as plt\n")
        v.write("import reaction_equations\n")
        v.write("if __name__ == \"__main__\":\n")
        #v.write("\ttempdata={}\n")
        v.write("\tmodel = kinetics.Model(logging=False)\n")
        #print(parameter_dict)
        #with open("reaction_mapping.txt","w") as mapping:
        
        for key in returned.keys():
            #if("medium" not in key):
            # if(key in tagged):
                # v.write("\t%s = reaction_equations.%s_Reaction(tempdata=tempdata)\n"%(key,key))
            # else:
            # if('vE' in key):
                # mapnameguess=key.replace("vE_","").upper()
                # print(mapnameguess)
                # mapping.write("%s\t%s\n"%(key,mapnameguess))
            if(key in mapping_names.keys()):
                v.write("\t%s = reaction_equations.%s_Reaction(metabolic_map_name=\"%s\")\n"%(key,key,mapping_names[key]))
            else:
                v.write("\t%s = reaction_equations.%s_Reaction()\n"%(key,key))
            v.write("\t%s.parameters={"%key)
            x=False
            for param in returned[key][1]:
                try:
                    value=parameter_dict[param]
                    v.write("'%s': %.10f,"%(param,value))
                    x=True
                    
                except:
                    x=False
            if(x):
                v.seek(v.tell()-1, os.SEEK_SET)
                v.truncate()
            v.write("}\n")
            v.write("\tmodel.append(%s)\n"%key)
        v.write("\tmodel.set_time(0, 300, 1000)\n\n")
        v.write("\tmodel.species = {")
        for key in init.keys():
            v.write("\"%s\" : %f,\n\t\t"%(key,init[key]))
        v.seek(v.tell()-5, os.SEEK_SET)
        v.truncate()
        v.write("\n\t\t}\n")
        v.write("\n\tmodel.setup_model()")
        v.write("\n\tmodel.run_model()")
        #v.write("\n\tmapping_names={")
        #for i in mapping_names.keys():
        #    v.write("\"%s\" : \"%s\",\n\t\t"%(i,mapping_names[i]))
        #v.seek(v.tell()-5, os.SEEK_SET)
        #v.truncate()
        #v.write("\n\t}")
            
if __name__ == "__main__":
    desired_reactions=["Glk","Pgi","Pfk","Fbp","Fba"]
    param_dict,unidict=parse_txt_file_parameters("kuratadata.txt")
    (returned,tagged)=reaction_class_creator_new("model_states.txt","reaction_equations.txt","substrate_list.txt","reaction_equations.py","model_variables.txt")
    init=get_initial_conditions("initial_conditions.txt")
    mapping_names=get_mapping_names("reaction_mapping.txt")
    write_running_code("run_metabolic_model.py", mapping_names, unidict,returned,init,tagged)