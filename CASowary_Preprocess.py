#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 15:04:03 2020

@author: xander
"""

import argparse, random, os
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
    target_Genes = {}
    trx_Seq = ''
    trx_Name = ''
    missed_Genes = []
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
    if (len(target_Genes) > 0):
        gene_Trx = Get_Gene_Transcripts()
        with open (output_File, 'w') as file:
            for i in gene_Name:
                if (i in gene_Trx.keys()):
                    for j in gene_Trx[i]:
                        for k in range(len(target_Genes[j]) - 27):
                            file.write('>' + i + '_' + str(k + 1) + '|' + str(j) + '\n')
                            file.write(target_Genes[j][k : k + 28] + '\n')
                else:
                    missed_Genes.append(i)
            file.close()
        print('Results saved as: ' + output_File + '\n')
        if (len(missed_Genes) > 0):
            for i in missed_Genes:
                print ('No results found for gene name: ' + i)
            print ('\n')
    else:
        print('No results available for given gene names.')

arg_Parser = argparse.ArgumentParser(description = 'Creates a fasta file containing all possible guides for all transcripts mapped to the input gene names.')
arg_Parser.add_argument('input', help = 'A gene name or list of gene names seperated by commas or a csv file containing gene names.')
arg_Parser.add_argument('-output', '-o', help = 'Name of the outputed fasta file name.')
fun_Args = arg_Parser.parse_args()


input_List = []
if (fun_Args.output == None):
    output_File = '{:010d}'.format(random.randint(0, 10000000000)) + '.fasta'
else:
    output_File = fun_Args.output + '.fasta'
if ('.csv' in fun_Args.input):
    with open (fun_Args.input) as file:
        text = file.readlines()
        file.close()
    for line in text:
        input_List += line.strip().split(',')
else:
    input_List += fun_Args.input.split(',')
for i in range(len(input_List)):
    input_List[i] = input_List[i].upper()
Determine_Guides(input_List)