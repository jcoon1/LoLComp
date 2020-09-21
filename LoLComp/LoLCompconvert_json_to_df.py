import pandas as pd
import json

def convert_json_to_df(trials_json):
    
    df = pd.DataFrame(trials_json['trials'])

    for key in trials_json['subject_information'].keys():
        df[key] = trials_json['subject_information'][key]

    for key in trials_json['system'].keys():
        df[key] = trials_json['system'][key]

    # for key in trials_json['subject_expertise'].keys():
    #   df[key] = trials_json['subject_expertise'][key]     

    df['time_in_minutes'] = trials_json['time_in_minutes']
    # df['workerID'] = trials_json['workerID']
    
    return df

#load data and name it json_trials
with open('LoLCompResults.json') as json_file:
    
    # add "[" as first element, remove last comma from end of string, then add "]" to end of string
    # json_file is an _io.textiowrapper object and the .read() method is converting it to a string
    json_file_string = '[' + json_file.read()[:-1] + ']'

    # json.loads expects a string as an iput 
    json_trials = json.loads(json_file_string)

df_trials = pd.DataFrame() #this creates an empty dataframe
participant_id = 1
for x in json_trials:
  df_ = convert_json_to_df(x)
  df_['participant_id'] = participant_id

  df_trials = pd.concat([df_trials, df_], sort=False)

  participant_id += 1

df_trials.to_csv("LoLCompResults.csv")

# print(json_trials)