import json
import os
import configparser

my_path = os.getcwd()

if __name__=='__main__':
    parser = configparser.ConfigParser()
    pars = {}
    parser.read(os.path.join(my_path, "params", "racecar1.ini"))
        
    pars["ggv_file"] = json.loads(parser.get('GENERAL_OPTIONS', 'ggv_file'))
    pars["ax_max_machines_file"] = json.loads(parser.get('GENERAL_OPTIONS', 'ax_max_machines_file'))
    pars["stepsize_opts"] = json.loads(parser.get('GENERAL_OPTIONS', 'stepsize_opts'))
    pars["reg_smooth_opts"] = json.loads(parser.get('GENERAL_OPTIONS', 'reg_smooth_opts'))
    pars["veh_params"] = json.loads(parser.get('GENERAL_OPTIONS', 'veh_params'))
    pars["vel_calc_opts"] = json.loads(parser.get('GENERAL_OPTIONS', 'vel_calc_opts'))
    pars["optim_opts"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'optim_opts_shortest_path'))
    pars["optim_opts"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'optim_opts_mincurv'))
    
    pars["curv_calc_opts"] = json.loads(parser.get('GENERAL_OPTIONS', 'curv_calc_opts'))
    pars["optim_opts"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'optim_opts_mintime'))
    pars["vehicle_params_mintime"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'vehicle_params_mintime'))
    pars["tire_params_mintime"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'tire_params_mintime'))
    pars["pwr_params_mintime"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'pwr_params_mintime'))
    
    # modification of mintime options/parameters
    pars["vehicle_params_mintime"]["wheelbase"] = (pars["vehicle_params_mintime"]["wheelbase_front"]
                                                       + pars["vehicle_params_mintime"]["wheelbase_rear"])
    
    for vals in pars:
        print(vals)
    
    #print(pars["veh_params"]["v_max"])
    try:
        choice1=input('Select Option to modify (str): ')
        pars[choice1]
        if (len(pars[choice1].split(','))==1) :
            print("Current Setting")
            print(choice1 +' = ' + pars[choice1])
        else:
            print("HERE?")
            print(pars[choice1].split(','))
            for val in pars[choice1].split(','):
                print (list)
                
    except:
        print("Just enter something valid here")
        
    return True
