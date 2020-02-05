#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 10:46:46 2018

"""

import pandas as pd
import numpy as np
import os
import pymatgen as mg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import itertools

ruta=input('Set CIF-files address: '+'\n')

if not ruta:
    ruta='./'
print('Address set to ',ruta,'\n')

def Positions(ruta=ruta, archivo='1010902', saveas='generated_file'):
    

    '''Se comienza por cargar el archivo y obtener las operaciones de simetria y la
    cristalografia del archivo'''
    
    sustancia=mg.Structure.from_file(str(ruta)+str(archivo)+'.cif')
    sga=SpacegroupAnalyzer(sustancia)
    archivo=sga.get_conventional_standard_structure()
    text=str(archivo)
    
    sitios=int(text.split('\n')[4].split('(')[1].split(')')[0])
            
    abc=[float(item) for item in list(filter(None,text.split('\n')[2].split(':')[1].split(' ')))]
    angles=[float(item) for item in list(filter(None,text.split('\n')[3].split(':')[1].split(' ')))]
    lista=text.split('\n')[-sitios:]
        
    newlist=[list(filter(None,line.split(' '))) for line in lista]
    '''
    for line in lista:
        newlist.append(list(filter(None,line.split(' '))))
    '''
    newlist=np.asarray(newlist)
        
    motif=pd.DataFrame(newlist)[[1,2,3,4]]
    motif[2]=[float(i) for i in motif[2].values]
    motif[3]=[float(i) for i in motif[3].values]
    motif[4]=[float(i) for i in motif[4].values]
    
    '''A continuacion se crea un archivo que contiene los parametros de
    red y los sitios en coordenadas fraccionarias'''
        
    motif.columns = np.arange(len(motif.columns))
    
    xyz=motif.iloc[:,1:4].values

    volumen=abc[0]*abc[1]*abc[2]*np.sqrt(1-(np.cos(np.deg2rad(angles[0])))**2-(np.cos(np.deg2rad(angles[1])))**2-(np.cos(np.deg2rad(angles[2])))**2+2*np.cos(np.deg2rad(angles[0]))*np.cos(np.deg2rad(angles[1]))*np.cos(np.deg2rad(angles[2])))

    #La variable matrix convierte las coordenadas relativas a coordenadas absolutas en un sistema cartesiano
    matrix=np.array([[abc[0],abc[1]*np.cos(np.deg2rad(angles[2])),abc[2]*np.cos(np.deg2rad(angles[1]))],
                      [0,abc[1]*np.sin(np.deg2rad(angles[2])),abc[2]*(np.cos(np.deg2rad(angles[0]))-np.cos(np.deg2rad(angles[1]))*np.cos(np.deg2rad(angles[2])))/np.sin(np.deg2rad(angles[2]))],
                      [0,0,volumen/(abc[0]*abc[1]*np.sin(np.deg2rad(angles[2])))]])
    
    positions=np.round(np.matmul(xyz,matrix),4)
    positions=pd.DataFrame(positions)
    positions=positions.rename(columns={0:'x',1:'y',2:'z'})
    positions['element']=motif[0]
    positions['element']=positions['element'].map(lambda x: x.lstrip('+-').rstrip('0123456789:.'))
    positions['element']=positions['element'].map(lambda x: x.lstrip('0123456789:.').rstrip('+-'))
    positions['element']=positions['element'].str.replace('\d+','')

    positions=positions[['element','x','y','z']]

    positions=positions.round(4)
    
    atoms=len(positions)

    positions.to_csv(str(saveas)+'.csv', header=None, index=None, sep='\t')
    
    with open(str(saveas)+'.xyz', 'w') as file:
        file.write(str(atoms)+'\n'+'\n')
        file.write(positions.to_string(header=None, index=None, col_space=0))
        file.close()
    
    return positions

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

def UnitCell(ruta=ruta, archivo='1010902', saveas='generated_file'):
    
    '''Se comienza por cargar el archivo y obtener las operaciones de simetria y la
    cristalografia del archivo'''
    
    sustancia=mg.Structure.from_file(str(ruta)+str(archivo)+'.cif')
    sga=SpacegroupAnalyzer(sustancia)
    archivo=sga.get_conventional_standard_structure()
    text=str(archivo)
    '''A continuacion se crea un archivo que contiene los parametros de
    red y los sitios en coordenadas fraccionarias'''
    
    sitios=int(text.split('\n')[4].split('(')[1].split(')')[0])
    
        
    abc=[float(item) for item in list(filter(None,text.split('\n')[2].split(':')[1].split(' ')))]
    angles=[float(item) for item in list(filter(None,text.split('\n')[3].split(':')[1].split(' ')))]
    lista=text.split('\n')[-sitios:]
    
    for i in range(len(lista)):
        print(lista[i])
    newlist=[list(filter(None,line.split(' '))) for line in lista]

    newlist=np.asarray(newlist)
    print(newlist)    
    motif=pd.DataFrame(newlist)[[1,2,3,4]]
    motif[2]=[float(i) for i in motif[2].values]
    motif[3]=[float(i) for i in motif[3].values]
    motif[4]=[float(i) for i in motif[4].values]
    
    motif.columns = np.arange(len(motif.columns))
    
    motif[4]=1*(motif[1].values == 0)+1*(motif[2].values == 0)+1*(motif[3].values == 0) #Se crea la columna de ayuda 'suma' para ver que atomos ocupan posiciones igual a cero
    
    motif_three=motif[motif[4] == 3].iloc[:,0:4].reset_index(drop=True) #Se elimina la columna de ayuda 'suma'
    motif_two=motif[motif[4] == 2].iloc[:,0:4].reset_index(drop=True) #Se elimina la columna de ayuda 'suma'
    motif_one=motif[motif[4] == 1].iloc[:,0:4].reset_index(drop=True) #Se elimina la columna de ayuda 'suma'
    motif=motif.iloc[:,0:4].reset_index(drop=True)

    #abc=np.loadtxt('abc.csv')
    #angles=np.loadtxt('angulos.csv')

    volumen=abc[0]*abc[1]*abc[2]*np.sqrt(1-(np.cos(np.deg2rad(angles[0])))**2-(np.cos(np.deg2rad(angles[1])))**2-(np.cos(np.deg2rad(angles[2])))**2+2*np.cos(np.deg2rad(angles[0]))*np.cos(np.deg2rad(angles[1]))*np.cos(np.deg2rad(angles[2])))

    #La variable matrix convierte las coordenadas relativas a coordenadas absolutas en un sistema cartesiano
    matrix=np.array([[abc[0],abc[1]*np.cos(np.deg2rad(angles[2])),abc[2]*np.cos(np.deg2rad(angles[1]))],
                      [0,abc[1]*np.sin(np.deg2rad(angles[2])),abc[2]*(np.cos(np.deg2rad(angles[0]))-np.cos(np.deg2rad(angles[1]))*np.cos(np.deg2rad(angles[2])))/np.sin(np.deg2rad(angles[2]))],
                      [0,0,volumen/(abc[0]*abc[1]*np.sin(np.deg2rad(angles[2])))]])

    if len(motif_one) != 0:
        
        positions_one=motif_one.replace(0,1).iloc[:,1:].values
        #positions_one=np.round(np.matmul(xyz_one,matrix),4)
        positions_one=pd.DataFrame(positions_one)
        positions_one=positions_one.rename(columns={0:1,1:2,2:3})
        positions_one[0]=motif_one[0]
        motif=motif.append(positions_one, ignore_index=True)
        
    tras0=np.array(list(itertools.product([0,1],repeat=3)))[1:7]
    tras1=np.array(list(itertools.product([0,1],repeat=3)))[1:]
    
        
    if len(motif_two) != 0:
        for vector in tras0:
            
            positions_two=np.add(motif_two.iloc[:,1:].values,vector)
            data=pd.DataFrame(positions_two)
            data=data.rename(columns={0:1,1:2,2:3})
            data[0]=motif_two[0]
            data=data[[0,1,2,3]]
            data=data[data <= 1].dropna()
            data=data.reset_index(drop=True)
            motif=motif.append(data, ignore_index=True)

    if len(motif_three) != 0:
        for vector in tras1:
            
            positions_three=np.add(motif_three.iloc[:,1:].values,vector)
            positions_three=pd.DataFrame(positions_three)
            data=pd.DataFrame(positions_three)
            data=data.rename(columns={0:1,1:2,2:3})
            data[0]=motif_three[0]
            data=data[[0,1,2,3]]
            data=data[data <= 1].dropna()
            data=data.reset_index(drop=True)
            motif=motif.append(data, ignore_index=True)
    
    positions=np.round(np.matmul(motif.iloc[:,1:].values,matrix),4)
    positions=pd.DataFrame(positions)
    positions=positions.rename(columns={0:'x',1:'y',2:'z'})
    positions['element']=motif[0]
    
    positions['element']=positions['element'].map(lambda x: x.lstrip('+-').rstrip('0123456789:.'))
    positions['element']=positions['element'].map(lambda x: x.lstrip('0123456789.:').rstrip('+-:.'))
    positions['element']=positions['element'].str.replace('\d+','')
    
    
    positions=positions[['element','x','y','z']]
    positions=positions.round(4)
    atoms=len(positions)
    
    positions.to_csv(str(saveas)+'.csv', header=None, index=None, sep='\t')
    
    with open(str(saveas)+'.xyz', 'w') as file:
        file.write(str(atoms)+'\n'+'\n')
        file.write(positions.to_string(header=None, index=None, col_space=0))
        file.close()
    
    return positions

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

def CrystalMaker(ruta=ruta,archivo='1010902', n=1, saveas='generated_file'):
    
    '''Se comienza por cargar el archivo y obtener las operaciones de simetria y la
    cristalografia del archivo'''
    
    sustancia=mg.Structure.from_file(str(ruta)+str(archivo)+'.cif')
    sga=SpacegroupAnalyzer(sustancia)
    archivo=sga.get_conventional_standard_structure()
    text=str(archivo)
    '''A continuacion se crea un archivo que contiene los parametros de
    red y los sitios en coordenadas fraccionarias'''
    
    sitios=int(text.split('\n')[4].split('(')[1].split(')')[0])
    
        
    abc=[float(item) for item in list(filter(None,text.split('\n')[2].split(':')[1].split(' ')))]
    angles=[float(item) for item in list(filter(None,text.split('\n')[3].split(':')[1].split(' ')))]
    lista=text.split('\n')[-sitios:]

    newlist=[list(filter(None,line.split(' '))) for line in lista]
    newlist=np.asarray(newlist)
            
    motif=pd.DataFrame(newlist)[[1,2,3,4]]
    motif[2]=[float(i) for i in motif[2].values]
    motif[3]=[float(i) for i in motif[3].values]
    motif[4]=[float(i) for i in motif[4].values]
    
    motif.columns = np.arange(len(motif.columns))
    
    motif[4]=1*(motif[1].values == 0)+1*(motif[2].values == 0)+1*(motif[3].values == 0) #Se crea la columna de ayuda 'suma' para ver que atomos ocupan posiciones igual a cero
    
    motif_three=motif[motif[4] == 3].iloc[:,0:4].reset_index(drop=True) #Se elimina la columna de ayuda 'suma'
    motif_two=motif[motif[4] == 2].iloc[:,0:4].reset_index(drop=True) #Se elimina la columna de ayuda 'suma'
    motif_one=motif[motif[4] == 1].iloc[:,0:4].reset_index(drop=True) #Se elimina la columna de ayuda 'suma'
    motif=motif.iloc[:,0:4].reset_index(drop=True)
    print(motif)
    volumen=abc[0]*abc[1]*abc[2]*np.sqrt(1-(np.cos(np.deg2rad(angles[0])))**2-(np.cos(np.deg2rad(angles[1])))**2-(np.cos(np.deg2rad(angles[2])))**2+2*np.cos(np.deg2rad(angles[0]))*np.cos(np.deg2rad(angles[1]))*np.cos(np.deg2rad(angles[2])))

    #La variable matrix convierte las coordenadas relativas a coordenadas absolutas en un sistema cartesiano
    matrix=np.array([[abc[0],abc[1]*np.cos(np.deg2rad(angles[2])),abc[2]*np.cos(np.deg2rad(angles[1]))],
                      [0,abc[1]*np.sin(np.deg2rad(angles[2])),abc[2]*(np.cos(np.deg2rad(angles[0]))-np.cos(np.deg2rad(angles[1]))*np.cos(np.deg2rad(angles[2])))/np.sin(np.deg2rad(angles[2]))],
                      [0,0,volumen/(abc[0]*abc[1]*np.sin(np.deg2rad(angles[2])))]])

    if len(motif_one) != 0:
        
        positions_one=motif_one.replace(0,1).iloc[:,1:].values
        #positions_one=np.round(np.matmul(xyz_one,matrix),4)
        positions_one=pd.DataFrame(positions_one)
        positions_one=positions_one.rename(columns={0:1,1:2,2:3})
        positions_one[0]=motif_one[0]
        motif=motif.append(positions_one, ignore_index=True)
        
    tras0=np.array(list(itertools.product([0,1],repeat=3)))[1:7]
    tras1=np.array(list(itertools.product([0,1],repeat=3)))[1:]
    
        
    if len(motif_two) != 0:
        for vector in tras0:
            
            positions_two=np.add(motif_two.iloc[:,1:].values,vector)
            data=pd.DataFrame(positions_two)
            data=data.rename(columns={0:1,1:2,2:3})
            data[0]=motif_two[0]
            data=data[[0,1,2,3]]
            data=data[data <= 1].dropna()
            data=data.reset_index(drop=True)
            motif=motif.append(data, ignore_index=True)

    if len(motif_three) != 0:
        for vector in tras1:
            
            positions_three=np.add(motif_three.iloc[:,1:].values,vector)
            positions_three=pd.DataFrame(positions_three)
            data=pd.DataFrame(positions_three)
            data=data.rename(columns={0:1,1:2,2:3})
            data[0]=motif_three[0]
            data=data[[0,1,2,3]]
            data=data[data <= 1].dropna()
            data=data.reset_index(drop=True)
            motif=motif.append(data, ignore_index=True)
    
    if n != 1:
        
        rows=len(motif)
        
        print('Each lattice parameter is increased by '+str(n)+' time(s)'+'\n')

        multiple=[m for m in range(n)]
        
        traslation=np.array(list(itertools.product(multiple, repeat=3)))[1:]
        
        for vector in traslation:
            data=np.add(motif.iloc[:rows,1:].values, vector)
            data=pd.DataFrame(data)
            data=data.rename(columns={0:1,1:2,2:3})
            data[0]=motif.iloc[:rows,0].values
            data=data[[0,1,2,3]]
            data=data.reset_index(drop=True)
            motif=motif.append(data, ignore_index=True)
                        
        motif=motif.round(4)
        data=data[data <= n].dropna()
        motif=motif.drop_duplicates() 
        motif=motif.reset_index(drop=True)
    print(motif)
    positions=np.round(np.matmul(motif.iloc[:,1:].values,matrix),4)
    positions=pd.DataFrame(positions)
    positions=positions.rename(columns={0:'x',1:'y',2:'z'})
    positions['element']=motif[0]
        
    positions['element']=positions['element'].map(lambda x: x.lstrip('+-').rstrip('0123456789'))
    positions['element']=positions['element'].map(lambda x: x.lstrip('0123456789').rstrip('+-'))
    positions['element']=positions['element'].str.replace('\d+','')

    positions=positions[['element','x','y','z']]
    positions=positions.round(4)

    atoms=len(positions)

    positions.to_csv(str(saveas)+'.csv', header=None, index=None, sep='\t')
    
    with open(str(saveas)+'.xyz', 'w') as file:
        file.write(str(atoms)+'\n'+'\n')
        file.write(positions.to_string(header=None, index=None, col_space=0))
        file.close()
    
    
    return positions

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'''
print('This is CrystalMaker, a program to generate xyz files.','\n',
      '\n','This program has three functions:','\n','\n',
      '1. Positions: This function specifies the atoms in the asymmetric cell','\n',
      '2. UnitCell: This function speciies the atoms in one unit cell'+'\n',
      '3. CrystalMaker: This function allows to grow the crystal in all dimentions by an integer factor n',
      '\n','\n','All positions are described in cartesian coordinates.','\n','\n',
      'After execution of a function, the program generates a csv - file and a xyz - file','\n','\n',
      'To leave the program, just type exit when this one asks you for a function.',
      'Please press Enter...','\n','\n')

order='start'

while order != 'exit':
    
    order = input('Which function would you like to use? [Positions/UnitCell/CrystalMaker]:'+'\n'+'\n')
    
    if order != 'exit':

        if order == 'Positions':
            try:
                archivo = input('Give me the cif number:'+'\n'+'\n')
                saveas = input('Give me a name to save the generated files:'+'\n'+'\n')
        
                Positions(archivo=archivo, saveas=saveas)
            except:
                print('There is an invalid entry. Please check typing!!!')    
            
        elif order == 'UnitCell':
            try:
                archivo = input('Give me the cif number:'+'\n'+'\n')
                saveas = input('Give me a name to save the generated files:'+'\n'+'\n')
        
                UnitCell(archivo=archivo, saveas=saveas)
            except:
                print('There is an invalid entry. Please check typing!!!')
        
        elif order == 'CrystalMaker':
            try:
                archivo = input('Give me the cif number:'+'\n'+'\n')
                saveas = input('Give me a name to save the generated files:'+'\n'+'\n')
                n = int(input('Multipy all dimensions by: '))
        
                CrystalMaker(archivo=archivo,n=n,saveas=saveas)
            except:
                print('There is an invalid entry. Please check typing!!!')
                
        else:
            print('Unknown function. Please check typing!!!'+'\n'+'\n')
            
    else:
        quit()
        
        

#'''    
