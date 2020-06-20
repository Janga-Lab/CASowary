#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 10:50:00 2019

@author: xander
"""

import numpy as np, os, sys
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

def Get_Transcript_Lengths():
    trx_Length = {}
    with open ('Gene_Information.txt') as file:
        text = file.readlines()
        for i in range(len(text)):
            line = text[i]
            if (line[0] == '>'):
                line = line.strip().split('|')
                trx_Length[line[1]] = int(line[7])
    file.close()
    return trx_Length

def Get_Guide_Information (fasta):
    guide_Name = []
    target_Trx = []
    guide_Sequences = []
    with open (fasta) as file:
        text = file.readlines()
        for i in range(len(text)):
            line = text[i].strip()
            if (line[0] == '>'):
                line = line.split('|')
                guide_Name.append(line[0][1:])
                target_Trx.append(line[1])
            else:
                sequence = line
                sequence = sequence[::-1]
                guide = ''
                for j in range(len(sequence)):
                    if (sequence[j] == 'A'):
                        guide += 'T'
                    if (sequence[j] == 'C'):
                        guide += 'G'
                    if (sequence[j] == 'G'):
                        guide += 'C'
                    if (sequence[j] == 'T'):
                        guide += 'A'
                guide_Sequences.append(guide)
    file.close()
    return guide_Name, target_Trx, guide_Sequences

def Get_Sequence_Features (guide_Sequences, features = 'Sig'):
    file_Names = {'Sig':'Significant_k-mers.txt','2':'Significant_k-mers_2Plus.txt','3':'Significant_k-mers_3Plus.txt','Gini':'Gini_k-mers.txt','DT':'Decision_Tree_k-mers.txt'}
    if (features in file_Names.keys()):
        with open (file_Names[features]) as file:
                text = file.readlines()
                sequence_Data = np.zeros([len(guide_Sequences),len(text) + 4])
                for i in range(len(text)):
                    line = text[i].strip().split('\t')
                    kmer = line[0]
                    index = int(line[1])
                    length = len(kmer)
                    for j in range(len(guide_Sequences)):
                        guide = guide_Sequences[j]
                        if (i == 0):
                            sequence_Data[j,-4] = guide.count('A')/len(guide)
                            sequence_Data[j,-3] = guide.count('T')/len(guide)
                            sequence_Data[j,-2] = guide.count('C')/len(guide)
                            sequence_Data[j,-1] = guide.count('G')/len(guide)
                        if (guide[index:index + length] == kmer):
                            sequence_Data[j,i] = 1
        file.close()
        return sequence_Data
    else:
        print ('Invalid feature list name. Try Sig, 2, 3, Gini, or DT.')
        return None

def Get_Availability_Features (fasta, bed):
    guide_Name, target_Trx, target_Seq = Get_Guide_Information (fasta)
    trx_Length = Get_Transcript_Lengths()
    target_Start = []
    target_End = []
    target_Location = []
    target_Hits = []
    hit_Id = []
    with open (hits_File) as file:
        text = file.readlines()
        for i in range (len(text)):
            line = text[i].strip().split('\t')
            hit_Id.append(line[3])
        for i in range(len(guide_Name)):
            if ((guide_Name[i] + '|' + target_Trx[i]) in hit_Id):
                hit_Count = hit_Id.count(guide_Name[i] + '|' + target_Trx[i])
                index = hit_Id.index(guide_Name[i] + '|' + target_Trx[i])
                for j in range(index,index + hit_Count):
                    guide_Hits = []
                    guide_Start = []
                    guide_End = []
                    line = text[j].strip().split('\t')
                    if (target_Trx[i] in line[0]):
                        guide_Start.append(int(line[1]))
                        guide_End.append(int(line[2]))
                        location = (int(line[1]) + int(line[2])) / 2
                        guide_Hits.append(location / trx_Length[target_Trx[i]])
                target_Hits.append(hit_Count)
                if (len(guide_Hits) > 0):
                    choice = np.random.randint(0,len(guide_Hits))
                    target_Location.append(guide_Hits[choice])
                    target_Start.append(guide_Start[choice])
                    target_End.append(guide_End[choice])
                else:
                    target_Start.append(0)
                    target_End.append(0)
                    target_Location.append(0)
            else:
                target_Hits.append(0)
                target_Start.append(0)
                target_End.append(0)
                target_Location.append(0)
    file.close()
    occ_Start, occ_End = Get_Protein_Occupancy(bed)
    target_Occupancy = []
    for i in range(len(target_Trx)):
        if (target_Trx[i] in occ_Start.keys()):
            for j in range(len(occ_Start[target_Trx[i]])):
                if (target_Start[i] > occ_Start[target_Trx[i]][j]) and (target_End[i] < occ_End[target_Trx[i]][j]):
                    target_Occupancy.append(1)
                    break
                elif (target_Start[i] < occ_Start[target_Trx[i]][j]) and (target_End[i] > occ_Start[target_Trx[i]][j]) and (target_End[i] < occ_End[target_Trx[i]][j]):
                    target_Occupancy.append((target_End[i] - occ_Start[target_Trx[i]][j]) / len(target_Seq[i]))
                    break
                elif (target_Start[i] > occ_Start[target_Trx[i]][j]) and (target_Start[i] < occ_End[target_Trx[i]][j]) and (target_End[i] > occ_End[target_Trx[i]][j]):
                    target_Occupancy.append((occ_End[target_Trx[i]][j] - target_Start[i]) / len(target_Seq[i]))
                    break
                if (len(target_Occupancy) < i + 1):
                    target_Occupancy.append(0)
        else:
            target_Occupancy.append(0)
    return np.array([target_Hits, target_Location, target_Occupancy]).T

def Sort_Guide_Hits(hits):
    bed_Data = []
    with open (hits, 'r') as file:
        text = file.readlines()
        for i in range(len(text)):
            line = text[i].strip().split('\t')
            bed_Data.append(line)
        bed_Data = np.array(bed_Data)
        bed_Data = bed_Data[bed_Data[:,3].argsort()]
        file.close()
    with open (hits,'w') as file:
        for i in range(bed_Data.shape[0]):
            file.write(bed_Data[i][0] + '\t' + bed_Data[i][1] + '\t' + bed_Data[i][2] + '\t' + bed_Data[i][3] + '\t' + bed_Data[i][4] + '\t' + bed_Data[i][5] + '\n')
        file.close()
    return text

def Get_All_Features (fasta, bed, features):
    guide_Information = Get_Guide_Information(fasta)
    sequence_Features = Get_Sequence_Features(guide_Information[2], features)
    availability_Features = Get_Availability_Features(fasta, bed)
    all_Features = np.append(sequence_Features, availability_Features, axis = 1)
    return all_Features

def Get_Protein_Occupancy (bed):
    occ_Start = {}
    occ_End = {}
    with open (bed) as file:
        text = file.readlines()
        for i in range(0,len(text)):
            line = text[i].strip().split('\t')
            if (len(line) == 10):
                start_Loc = i + 1
                break
        for i in range(start_Loc, len(text)):
            line = text[i].strip().split('\t')
            transcript = line[0].split('.')[0]
            if (transcript in occ_Start.keys()):
                occ_Start[transcript].append(int(line[1]) - 1)
                occ_End[transcript].append(int(line[2]) - 1)
            else:
                occ_Start[transcript] = [int(line[1]) - 1]
                occ_End[transcript] = [int(line[2]) - 1]
    file.close()
    return occ_Start, occ_End

def Get_Training_Data (features = 'Sig'):
    file_Names = {'Sig':'Training Data Sig.txt','2':'Training Data 2.txt','3':'Training Data 3.txt','Gini':'Training Data Gini.txt','DT':'Training Data Decision Tree.txt'}
    raw_Data = []
    if (features in file_Names.keys()):
        with open (file_Names[features]) as file:
            text = file.readlines()
            for i in range(1,len(text)):
                line = text[i].strip().split('\t')
                raw_Data.append(line)
            file.close()
            raw_Data = np.array(raw_Data)
            raw_Data = raw_Data.astype('float')
            train_Data = raw_Data[:,0:-1]
            train_Label = raw_Data[:,-1]
            return train_Data, train_Label
    else:
        print ('Invalid feature list name. Try Sig, 2, 3, Gini, or DT.')
        return None

def Create_Model(model, features = 'Sig'):
    model_Types = ['Forest','KNN','Linear','Poly','Sigmoid','Tree']
    if (model in model_Types):
        if (model == 'Forest'):
            from sklearn import ensemble
            model_Parameters = {'Sig':40,'2':45,'3':45,'Gini':45,'DT':40}
            clf = ensemble.RandomForestClassifier(n_estimators = model_Parameters[features])
        elif (model == 'KNN'):
            from sklearn import neighbors
            model_Parameters = {'Sig':3,'2':10,'3':8,'Gini':3,'DT':3}
            clf = neighbors.KNeighborsClassifier(n_neighbors = model_Parameters[features], weights = 'uniform')
        elif (model == 'Tree'):
            from sklearn import tree
            clf = tree.DecisionTreeClassifier()
        else:
            from sklearn import svm
            clf = svm.SVC(max_iter = 10000, kernel = model, gamma = 'scale')
        return clf
    else:
        print ('Invalid model type. Try Forest, KNN, Linear, Polyn, Sigmoid, or Tree.')
        return None

def Parse_Results (result_File):
    guide_Results = {}
    with open (result_File) as file:
        text = file.readlines()
        for i in range(1,len(text)):
            line = text[i].strip().split('\t')
            gene = line[0].split('_')[0]
            transcript = line[1]
            predict = int(line[3][0])
            confidence = float(line[4 + predict])
            if (gene not in guide_Results.keys()):
                guide_Results[gene] = {}
            if (transcript not in guide_Results[gene].keys()):
                guide_Results[gene][transcript] = {}
            if (predict not in guide_Results[gene][transcript].keys()):
                guide_Results[gene][transcript][predict] = {}
            if (confidence not in guide_Results[gene][transcript][predict].keys()):
                guide_Results[gene][transcript][predict][confidence] = []
            guide_Results[gene][transcript][predict][confidence].append(i)
    file.close()
    path = os.getcwd()
    for i in guide_Results.keys():
        os.chdir(path)
        os.mkdir('Results_' + i)
        for j in guide_Results[i].keys():
            os.chdir(path + '/Results_' + i)
            os.mkdir(j)
            os.chdir(path + '/Results_' + i + '/' + j)
            for k in guide_Results[i][j].keys():
                list_Keys = list(guide_Results[i][j][k].keys())
                list_Keys.sort(reverse = True)
                if (k == 0) or (k == 1):
                    with open ('Best_Guides.txt','w') as file:
                        if (k == 0):
                            file.write('High Quality Guides (Class 0)' + '\n')
                        else:
                            file.write('Quality Guides (Class 1)' + '\n')
                        for l in list_Keys:
                            for m in guide_Results[i][j][k][l]:
                                line = text[m].strip().split('\t')
                                file.write(line[0] + '\t' + '5`-'+line[2] + '-3`' + '\n')
                            file.write('\n')
                        file.close()
                elif (k == 3):
                    with open ('Worst_Guides.txt','w') as file:
                        for l in list_Keys:
                            for m in guide_Results[i][j][k][l]:
                                line = text[m].strip().split('\t')
                                file.write(line[0] + '\t' + '5`-'+line[2] + '-3`' + '\n')
                            file.write('\n')
                        file.close()
    os.chdir(path)


fasta_File = 'Input_Guides.fasta'
occupy_File = 'SRR1033461_transcriptome_peaks.xls'
hits_File = 'Target_Hits.sorted.bed'
guide_Name, target_Trx, target_Seq = Get_Guide_Information (fasta_File)
test_Data = Get_All_Features(fasta_File, occupy_File, 'Gini')
train_Data, train_Label = Get_Training_Data('Gini')
model = Create_Model('Tree')
model.fit(train_Data, train_Label)
results = model.predict(test_Data)
probability = model.predict_proba(test_Data)
predict_File = 'Target_Predictions.txt'

with open (predict_File, 'w') as file:
    file.write('Guide Name' + '\t' + 'Transcript Name' + '\t' + 'Guide Sequence' + '\t' + 'Location' + '\t' + 'Model Prediction' + '\t' + 'Class 0 Probability' + '\t' + 'Class 1 Probability' + '\t' + 'Class 2 Probability' + '\t' + 'Class 3 Probability' + '\n')
    for i in range (len(guide_Name)):
        file.write(guide_Name[i] + '\t' + target_Trx[i] + '\t' + target_Seq[i] + '\t' + str(test_Data[i][110]) + '\t' + str(results[i]) + '\t' + str(probability[i,0]) + '\t' + str(probability[i,1]) + '\t' + str(probability[i,2]) + '\t' + str(probability[i,3]) + '\n')
    file.close()
