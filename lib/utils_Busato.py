#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 12:15:06 2023

@author: ben
"""


sec_h=3600
sec_g=24*3600

def atmbc_with_suspension(HSPATM,IETO,
                          zero,
                          inputs_atmbc,
                          nnod3d,nn,nn2d,
                          areanod2d,
                          x,y,z,
                          time,
                          ):
    
    counter=0
    totarea=0
    
    # Write atmbc file
    with open('atmbc', 'w') as file:
        # Count nodes for atmbc
        
        file.write(f"{HSPATM} {IETO} HSPATM  IETO\n")
        file.write(f"{zero:.2f} TIME\n")
        
        
        if inputs_atmbc['giorni_stop'] != 0 and inputs_atmbc['giorni_stop'] != 0:
            for i in range(nnod3d):
                if z[i] == zero:
                    if inputs_atmbc['xmin'] < x[i] < inputs_atmbc['xmax'] and inputs_atmbc['ymin'] < y[i] < inputs_atmbc['ymax']:
                        counter += 1
                        # print(counter)
                        if nn[i] == nn2d[i]:
                            print(counter)
                            # file.write(f"{nn[i]} {nn2d[i]} {areanod2d[i]}\n")
                            totarea += areanod2d[i]
            
            for i in range(nnod3d):
                if z[i] == zero:
                    file.write(f"{zero:.2f}\n")
    
    
            for k in range(inputs_atmbc['giorni_irrigazione']):
                file.write(f"{time:.2f} TIME\n")
                for i in range(nnod3d):
                    if z[i] == zero:
                        if inputs_atmbc['xmin'] < x[i] < inputs_atmbc['xmax'] and inputs_atmbc['ymin'] < y[i] < inputs_atmbc['ymax']:
                            file.write(f"{inputs_atmbc['flux'] / totarea}\n")
                        else:
                            file.write(f"{zero:.2f}\n")
    
                t = time + 5 * sec_h + 1
                file.write(f"\n{t:.2f} TIME\n")
                for i in range(nnod3d):
                    if z[i] == zero:
                        file.write(f"{zero:.2f}\n")
                file.write("\n")
                time += sec_g
                t += sec_g
        
        
            time += sec_g * inputs_atmbc['giorni_stop']
            file.write(f"{time:.2f} TIME\n")
            for i in range(nnod3d):
                if z[i] == zero:
                    if inputs_atmbc['xmin'] < x[i] < inputs_atmbc['xmax'] and inputs_atmbc['ymin'] < y[i] < inputs_atmbc['ymax']:
                        file.write(f"{inputs_atmbc['flux'] / totarea}\n")
                    else:
                        file.write(f"{zero:.2f}\n")
    
            t += sec_g * inputs_atmbc['giorni_stop']
            file.write(f"\n{t:.2f} TIME\n")
            for i in range(nnod3d):
                if z[i] == zero:
                    file.write(f"{zero:.2f}\n")
            file.write("\n")
            time += sec_g
            t += sec_g
            
            for k in range(inputs_atmbc['giorni_irrigazione_2'] - 1):
                file.write(f"{time:.2f} TIME\n")
                for i in range(nnod3d):
                    if z[i] == zero:
                        if inputs_atmbc['xmin'] < x[i] < inputs_atmbc['xmax'] and inputs_atmbc['ymin'] < y[i] < inputs_atmbc['ymax']:
                            file.write(f"{inputs_atmbc['flux'] / totarea}\n")
                        else:
                            file.write(f"{zero:.2f}\n")
            
                t = time + 5 * sec_h + 1
                file.write(f"\n{t:.2f} TIME\n")
                for i in range(nnod3d):
                    if z[i] == zero:
                        file.write(f"{zero:.2f}\n")
                file.write("\n")
                time += sec_g
                t += sec_g
    