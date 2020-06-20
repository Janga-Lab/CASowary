#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 15:04:03 2020

@author: xander
"""

import os, sys
os.path.join(os.path.dirname(__file__))

def Get_Gene_Transcripts():
    gene_Trx = {}
    with open ('Gene_Information.txt') as file:
        text = file.readlines()
        for i in range(len(text)):
            line = text[i]
            if (line[0] == '>'):
                line = line.strip().split('|')
                if (line[2] in gene_Trx.keys()):
                    gene_Trx[line[2]].append(line[1])
                else:
                    gene_Trx[line[2]] = [line[1]]
    file.close()
    return gene_Trx

def Determine_Guides (gene_Name):
    if (type(gene_Name) == str):
        gene_Name = list(gene_Name)
    elif (type(gene_Name) != list):
        print ('Invalid input type. Please run using string of gene name, or list of gene names.')
    else:
        target_Genes = {}
        trx_Seq = ''
        trx_Name = ''
        with open ('Gene_Information.txt') as file:
            text = file.readlines()
            for i in range(len(text)):
                line = text[i]
                if (line[0] == '>'):
                    if (trx_Seq != ''):
                        target_Genes[trx_Name] = trx_Seq
                        trx_Seq = ''
                    line = line.strip().split('|')
                    trx_Name = line[1]
                    trx_Gene = line[2]
                else:
                    if (trx_Gene in gene_Name):
                        line = line.strip()
                        trx_Seq += line
            target_Genes[trx_Name] = trx_Seq
        file.close()
        if (len(target_Genes) >0):
            gene_Trx = Get_Gene_Transcripts()
            with open ('Input_Guides.fasta', 'w') as file:
                for i in gene_Name:
                    for j in gene_Trx[i]:
                        for k in range(len(target_Genes[j]) - 27):
                            file.write('>' + i + '_' + str(k + 1) + '|' + str(j) + '\n')
                            file.write(target_Genes[j][k : k + 28] + '\n')
            file.close()
            print('Results saved as: \'Input_Guides.fasta\'\n')
        else:
            print('No results available for given gene names.')

input_File = sys.argv[1]
with open (input_File) as file:
    text = file.readlines()
    file.close()
input_List = []
for line in text:
    input_List += line.strip().split(',')
Determine_Guides(input_List)