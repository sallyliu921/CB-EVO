# Copyright (c) 2021, Programmable digital systems group, University of Toronto
# All rights reserved.

# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree. 


import re
import numpy as np
import datetime
from multiprocessing import Process, Manager
from subprocess import check_output
from collections import defaultdict

design_file = "../DRiLLS/benchmarks/arithmetic/multiplier.blif"
ccirc_binary = "Cgen/ccirc/ccirc"
# yosys_binary = 
# abc_binary = 


def shape_classifier(datas):
    
    gradients=np.diff(datas)
    maxima_num=0
    max_locations=[]
    count=0
    Threshold = 2* np.mean(datas)
    print("Threshold is " , Threshold)
    for i in gradients[:-1]:
        count+=1

        if ((i> Threshold/2 ) & (gradients[count]< Threshold/2 ) & (i != gradients[count]) ):
            maxima_num+=1
            max_locations.append(count)
    turning_points = {'maxima_number':maxima_num,'maxima_locations':max_locations}
    print (turning_points)
    exit()

    return shape_type


def ccirc_stats(design_file, ccirc_binary, stats):
    stats_file = '../../netlist.stats'
    #ccirc_command = design_file + " --partitions 1" + " --out " + stats_file

    
    try:
        proc = check_output([ccirc_binary, design_file, "--partitions" ," 1","--out","netlist.stats" ])
        readfile = "Cgen/ccirc/" + stats_file
        with open(readfile) as f:
            lines = f.readlines()   
        for line in lines:
            #Basic characteristic of Circuit (Feature count: 5)
            if 'Number_of_Nodes' in line:
                stats['Number_of_Nodes'] = int(line.strip().split()[-1])
            if 'Number_of_Edges' in line:
                stats['Number_of_Edges'] = int(line.strip().split()[-1])
            if 'Maximum_Delay' in line:
                stats['Maximum_Delay'] = float(line.strip().split()[-1])
            if 'Number_of_Combinational_Nodes' in line:
                stats['Number_of_Combinational_Nodes'] = float(line.strip().split()[-1])
            if ('Number_of_DFF' in line) and ('Number_of_DFF' not in stats):
                #Not useful in combination circuits, will always be zero
                stats['Number_of_DFF'] = float(line.strip().split()[-1])
            
            #Characteristic of Delay structure (Feature count: 7 * ?)
            if 'Node_shape' in line:
                stats['Node_shape'] = shape_classifier(np.array(line.strip().split()[2:-1],dtype=float))
            if 'Input_shape:' in line:
                stats['Input_shape:'] = np.array(line.strip().split()[2:-1])
            if 'Output_shape' in line:
                stats['Output_shape'] = np.array(line.strip().split()[2:-1])
            if 'Latched_shape' in line:
                #Not useful in combination circuits, will always be zero
                stats['Latched_shape'] = np.array(line.strip().split()[2:-1])
            if 'POshape' in line:
                stats['POshape'] = np.array(line.strip().split()[2:-1])
            if 'Edge_length_distribution' in line:
                # The same as Intra_cluster_edge_length_distribution for cluster count = 1
                stats['Edge_length_distribution'] = np.array(line.strip().split()[2:-1])
            if 'Fanout_distribution' in line:
                #Should use a different way to analyze the graph -- to be implemented
                stats['Fanout_distribution'] = np.array(line.strip().split()[2:-1])
            
            # #Characteristic of connections, (Feature count: 1 * ?)
            # if 'Intra_cluster_edge_length_distribution' in line:
            #     stats['Intra_cluster_edge_length_distribution'] = np.array(line.strip().split()[2:-1])

            #Characteristic of Fanout from Nodes (Feature Count: 14)
            if 'Maximum_fanout' in line:
                stats['Maximum_fanout'] = float(line.strip().split()[-1])
            if 'Number_of_high_degree_comb' in line:
                stats['Number_of_high_degree_comb'] = float(line.strip().split()[-1])
            if 'Number_of_high_degree_pi' in line:
                stats['Number_of_high_degree_pi'] = float(line.strip().split()[-1])
            if 'Number_of_high_degree_dff' in line:
                stats['Number_of_high_degree_dff'] = float(line.strip().split()[-1])
            if 'Number_of_10plus_degree_comb' in line:
                stats['Number_of_10plus_degree_comb'] = float(line.strip().split()[-1])
            if 'Number_of_10plus_degree_pi' in line:
                stats['Number_of_10plus_degree_pi'] = float(line.strip().split()[-1])
            if 'Avg_fanin' in line:
                stats['Avg_fanin' ] = float(line.strip().split()[-2])
                s = line.strip().split()[-1]
                stats['Std_fanin' ] = float(s[s.find("(")+1:s.find(")")])
            if 'Avg_fanin_comb' in line:
                stats['Avg_fanin_comb' ] = float(line.strip().split()[-2])
                s = line.strip().split()[-1]
                stats['Std_fanin_comb' ] = float(s[s.find("(")+1:s.find(")")])
            if 'Avg_fanin_pi' in line:
                stats['Avg_fanin_pi' ] = float(line.strip().split()[-2])
                s = line.strip().split()[-1]
                stats['Std_fanin_pi' ] = float(s[s.find("(")+1:s.find(")")])
            if 'Avg_fanin_dff' in line:
                stats['Avg_fanin_dff' ] = float(line.strip().split()[-2])
                s = line.strip().split()[-1]
                stats['Std_fanin_dff' ] = float(s[s.find("(")+1:s.find(")")])
            
    except Exception as e:
        print(e)
        return None
    return stats



def yosys_stats(design_file, yosys_binary, stats):
    yosys_command = "read_verilog " + design_file + "; stat"
    try:
        proc = check_output([yosys_binary, '-QT', '-p', yosys_command])
        lines = proc.decode("utf-8").split('\n')
        for line in lines:
            if 'Number of wires' in line:
                stats['number_of_wires'] = int(line.strip().split()[-1])
            if 'Number of public wires' in line:
                stats['number_of_public_wires'] = int(line.strip().split()[-1])
            if 'Number of cells' in line:
                stats['number_of_cells'] = float(line.strip().split()[-1])
            if '$and' in line:
                stats['ands'] = float(line.strip().split()[-1])
            if '$or' in line:
                stats['ors'] = float(line.strip().split()[-1])
            if '$not' in line:
                stats['nots'] = float(line.strip().split()[-1])
        
        # catch some design special cases
        if 'ands' not in stats:
            stats['ands'] = 0.0
        if 'ors' not in stats:
            stats['ors'] = 0.0
        if 'nots' not in stats:
            stats['nots'] = 0.0
    except Exception as e:
        print(e)
        return None
    return stats

def abc_stats(design_file, abc_binary, stats):    
    abc_command = "read_verilog " + design_file + "; print_stats"
    try:
        proc = check_output([abc_binary, '-c', abc_command])
        lines = proc.decode("utf-8").split('\n')
        for line in lines:
            if 'i/o' in line:
                ob = re.search(r'i/o *= *[0-9]+ */ *[0-9]+', line)
                stats['input_pins'] = int(ob.group().split('=')[1].strip().split('/')[0].strip())
                stats['output_pins'] = int(ob.group().split('=')[1].strip().split('/')[1].strip())
        
                ob = re.search(r'edge *= *[0-9]+', line)
                stats['edges'] = int(ob.group().split('=')[1].strip())

                ob = re.search(r'lev *= *[0-9]+', line)
                stats['levels'] = int(ob.group().split('=')[1].strip())

                ob = re.search(r'lat *= *[0-9]+', line)
                stats['latches'] = int(ob.group().split('=')[1].strip())
    except Exception as e:
        print(e)
        return None
    
    return stats

def extract_features(design_file, yosys_binary='yosys', abc_binary='abc'):
    '''
    Returns features of a given circuit as a tuple.
    Features are listed below
    '''
    manager = Manager()
    stats = manager.dict()
    p1 = Process(target=yosys_stats, args=(design_file, yosys_binary, stats))
    p2 = Process(target=abc_stats, args=(design_file, abc_binary, stats))
    p1.start()
    p2.start()
    p1.join()
    p2.join()

    p3 = Process(target=abc_stats, args=(design_file, abc_binary, stats))
    p3.start()
    p3.join()

    print(stats)


    # normalized features
    features = defaultdict(float)    
    
    

    return stats


if __name__ == '__main__':
    stats = defaultdict(list)    
    ccirc_stats(design_file, ccirc_binary, stats)
    print(stats)
