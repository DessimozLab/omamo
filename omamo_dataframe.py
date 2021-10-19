#!/usr/bin/python3

from sys import argv
import pandas as pd
import os, sys

def create_dataframe(directory_name):

	files = os.listdir(directory)
	csv_files = []
	for file in files:
    		if file.split(".")[-1]=="csv" and file.split(".")[0][-1] =="2":
            		csv_files.append(file)
	print('THERE ARE {} CSV FILES IN THE DIRECTORY'.format(len(csv_files)))
	df_list = []
	for file in csv_files:
    		df_list.append(pd.read_csv(directory+'/{}'.format(file), sep = '\t'))
	df = pd.concat(df_list, ignore_index = True)
	
	df = df.sort_values(by = ['GO ID', 'Score'],ignore_index = True, ascending = [True, False])
	
	del df['Unnamed: 0']
	
	return df


if __name__ == "__main__":
	
	'''ARG1 - address of the directory containing csv files'''
	
	directory = argv[1]

	df = create_dataframe(directory)

	df.to_csv('omamo_output_df.csv', sep = '\t')
	
	
